from flask import Blueprint, request, jsonify, current_app
from werkzeug.utils import secure_filename
# import os
import uuid
# import json   
import threading
from datetime import datetime
# from Bio import SeqIO
# from io import StringIO
# import tempfile
# import asyncio
from concurrent.futures import ThreadPoolExecutor
from src.visualizations import VisualizationEngine

from database.db_manager import DatabaseManager
from src.primer_designer import PrimerDesigner
from src.utils import allowed_file, parse_fasta_content, validate_parameters
from config import Config

api_bp = Blueprint('api', __name__)
db_manager = DatabaseManager()
executor = ThreadPoolExecutor(max_workers=Config.MAX_CONCURRENT_JOBS)

# Job status tracking
job_status = {}
job_lock = threading.Lock()

@api_bp.route('/design', methods=['POST'])
def design_primers():
    """Main endpoint for primer design requests"""
    try:
        # Generate unique job ID
        job_id = str(uuid.uuid4())
        
        # Parse request data
        if request.is_json:
            # JSON request with sequence data
            data = request.get_json()
            sequence = data.get('sequence', '')
            parameters = data.get('parameters', {})
            
            if not sequence:
                return jsonify({
                    'error': 'No sequence provided'
                }), 400
                
            # Format sequence data
            sequences = [{
                'name': 'target_sequence',
                'sequence': sequence.upper().replace(' ', '').replace('\n', '')
            }]
            job_type = 'json_input'
        else:
            # Form data with file upload
            sequences = []
            parameters = {}
            
            # Extract parameters from form
            for key, value in request.form.items():
                if key.startswith('param_'):
                    param_name = key.replace('param_', '')
                    try:
                        parameters[param_name] = float(value) if '.' in value else int(value)
                    except ValueError:
                        parameters[param_name] = value
            
            # Handle file upload
            if 'file' in request.files:
                file = request.files['file']
                if file.filename != '' and allowed_file(file.filename):
                    try:
                        # Read and parse FASTA content
                        content = file.read().decode('utf-8')
                        sequences = parse_fasta_content(content)
                        job_type = 'file_upload'
                    except Exception as e:
                        return jsonify({
                            'error': f'Error reading file: {str(e)}'
                        }), 400
                else:
                    return jsonify({
                        'error': 'Invalid file format. Please upload a FASTA file.'
                    }), 400
            else:
                return jsonify({
                    'error': 'No sequences provided'
                }), 400
        
        # Validate input
        if not sequences:
            return jsonify({
                'error': 'No valid sequences found'
            }), 400
        
        # Validate and set default parameters
        parameters = validate_parameters(parameters)
        
        # Create job in database
        input_data = {'sequences': sequences}
        db_manager.create_job(job_id, job_type, input_data, parameters)
        
        # Submit job for processing
        future = executor.submit(process_primer_design, job_id, sequences, parameters)
        
        # Track job status
        with job_lock:
            job_status[job_id] = {
                'status': 'pending',
                'created_at': datetime.now().isoformat(),
                'future': future
            }
        
        return jsonify({
            'job_id': job_id,
            'status': 'pending',
            'message': 'Primer design job submitted successfully',
            'estimated_time': len(sequences) * 10  # Rough estimate in seconds
        }), 202
        
    except Exception as e:
        current_app.logger.error(f"Error in primer design: {str(e)}")
        return jsonify({
            'error': 'Internal server error occurred'
        }), 500

@api_bp.route('/status/<job_id>', methods=['GET'])
def get_job_status(job_id):
    """Get job status and results"""
    try:
        # Check local job status first
        with job_lock:
            if job_id in job_status:
                local_status = job_status[job_id]
                if local_status['future'].done():
                    # Job completed, clean up
                    del job_status[job_id]
        
        # Get status from database
        job_data = db_manager.get_job(job_id)
        
        if not job_data:
            return jsonify({
                'error': 'Job not found'
            }), 404
        
        response_data = {
            'job_id': job_id,
            'status': job_data['status'],
            'created_at': job_data['created_at'],
            'updated_at': job_data['updated_at']
        }
        
        if job_data['status'] == 'completed':
            response_data['results'] = job_data['results']
            response_data['completed_at'] = job_data['completed_at']
        elif job_data['status'] == 'failed':
            response_data['error_message'] = job_data['error_message']
        
        return jsonify(response_data)
        
    except Exception as e:
        current_app.logger.error(f"Error getting job status: {str(e)}")
        return jsonify({
            'error': 'Error retrieving job status'
        }), 500

@api_bp.route('/batch', methods=['POST'])
def batch_process():
    """Handle batch processing of multiple files"""
    try:
        if 'files' not in request.files:
            return jsonify({
                'error': 'No files provided'
            }), 400
        
        files = request.files.getlist('files')
        if not files or all(f.filename == '' for f in files):
            return jsonify({
                'error': 'No files selected'
            }), 400
        
        # Validate file count
        if len(files) > 10:  # Limit batch size
            return jsonify({
                'error': 'Too many files. Maximum 10 files per batch.'
            }), 400
        
        # Extract parameters
        parameters = {}
        for key, value in request.form.items():
            if key.startswith('param_'):
                param_name = key.replace('param_', '')
                try:
                    parameters[param_name] = float(value) if '.' in value else int(value)
                except ValueError:
                    parameters[param_name] = value
        
        parameters = validate_parameters(parameters)
        
        # Process each file
        batch_id = str(uuid.uuid4())
        job_ids = []
        
        for i, file in enumerate(files):
            if file.filename != '' and allowed_file(file.filename):
                try:
                    # Parse file content
                    content = file.read().decode('utf-8')
                    sequences = parse_fasta_content(content)
                    
                    if sequences:
                        # Create individual job
                        job_id = f"{batch_id}_{i}"
                        input_data = {
                            'sequences': sequences,
                            'filename': secure_filename(file.filename),
                            'batch_id': batch_id
                        }
                        
                        db_manager.create_job(job_id, 'batch_processing', input_data, parameters)
                        
                        # Submit for processing
                        future = executor.submit(process_primer_design, job_id, sequences, parameters)
                        
                        with job_lock:
                            job_status[job_id] = {
                                'status': 'pending',
                                'created_at': datetime.now().isoformat(),
                                'future': future,
                                'batch_id': batch_id
                            }
                        
                        job_ids.append(job_id)
                        
                except Exception as e:
                    current_app.logger.error(f"Error processing file {file.filename}: {str(e)}")
        
        if not job_ids:
            return jsonify({
                'error': 'No valid files could be processed'
            }), 400
        
        return jsonify({
            'batch_id': batch_id,
            'job_ids': job_ids,
            'status': 'pending',
            'message': f'Batch processing started for {len(job_ids)} files'
        }), 202
        
    except Exception as e:
        current_app.logger.error(f"Error in batch processing: {str(e)}")
        return jsonify({
            'error': 'Internal server error occurred'
        }), 500

@api_bp.route('/batch/status/<batch_id>', methods=['GET'])
def get_batch_status(batch_id):
    """Get status of all jobs in a batch"""
    try:
        # Get all jobs for this batch
        conn = db_manager.db_path
        import sqlite3
        
        db_conn = sqlite3.connect(conn)
        cursor = db_conn.cursor()
        
        cursor.execute('''
            SELECT id, status, created_at, completed_at, error_message
            FROM jobs 
            WHERE id LIKE ? 
            ORDER BY id
        ''', (f"{batch_id}_%",))
        
        rows = cursor.fetchall()
        db_conn.close()
        
        if not rows:
            return jsonify({
                'error': 'Batch not found'
            }), 404
        
        jobs = []
        completed_count = 0
        failed_count = 0
        
        for row in rows:
            job_data = {
                'job_id': row[0],
                'status': row[1],
                'created_at': row[2],
                'completed_at': row[3],
                'error_message': row[4]
            }
            jobs.append(job_data)
            
            if row[1] == 'completed':
                completed_count += 1
            elif row[1] == 'failed':
                failed_count += 1
        
        overall_status = 'pending'
        if completed_count + failed_count == len(jobs):
            overall_status = 'completed'
        elif failed_count > 0:
            overall_status = 'partial_failure'
        
        return jsonify({
            'batch_id': batch_id,
            'overall_status': overall_status,
            'total_jobs': len(jobs),
            'completed': completed_count,
            'failed': failed_count,
            'pending': len(jobs) - completed_count - failed_count,
            'jobs': jobs
        })
        
    except Exception as e:
        current_app.logger.error(f"Error getting batch status: {str(e)}")
        return jsonify({
            'error': 'Error retrieving batch status'
        }), 500

@api_bp.route('/history', methods=['GET'])
def get_job_history():
    """Get job history"""
    try:
        limit = request.args.get('limit', 50, type=int)
        history = db_manager.get_job_history(limit)
        
        return jsonify({
            'jobs': history,
            'total': len(history)
        })
        
    except Exception as e:
        current_app.logger.error(f"Error getting job history: {str(e)}")
        return jsonify({
            'error': 'Error retrieving job history'
        }), 500

def process_primer_design(job_id, sequences, parameters):
    from app import app  # Import here to avoid circular import
    with app.app_context():
        print(f"[PRINT-DEBUG] Entered process_primer_design for job {job_id}")
        try:
            # Update status to processing
            db_manager.update_job_status(job_id, 'processing')
            
            # Log parameters used for this job
            current_app.logger.info(f"[DEBUG] Job {job_id} - Parameters used: {parameters}")
            
            # Initialize primer designer with parameters
            designer = PrimerDesigner()
            
            # Process each sequence
            all_results = []
            for seq_data in sequences:
                try:
                    # Design primers for this sequence
                    result = designer.design_primers(
                        sequence=seq_data['sequence'],
                        custom_params=parameters,
                        num_return=5
                    )
                    
                    if result.get('success', False) and result.get('primer_pairs'):
                        # Format results for storage
                        for pair in result['primer_pairs']:
                            # Log quality scores for this pair
                            left_q = pair['left_primer']['analysis']['quality_scores']
                            right_q = pair['right_primer']['analysis']['quality_scores']
                            overall_q = pair['overall_quality']
                            current_app.logger.info(
                                f"[DEBUG] Job {job_id} - Pair {pair['pair_index']} - "
                                f"Left quality: {left_q}, Right quality: {right_q}, Overall: {overall_q}"
                            )
                            primer_data = {
                                'sequence_name': seq_data.get('name', 'unknown'),
                                'forward_primer': pair['left_primer']['sequence'],
                                'reverse_primer': pair['right_primer']['sequence'],
                                'forward_tm': pair['left_primer']['analysis']['melting_temperature'],
                                'reverse_tm': pair['right_primer']['analysis']['melting_temperature'],
                                'forward_position': f"{pair['left_primer']['start']}-{pair['left_primer']['start'] + pair['left_primer']['length']}",
                                'reverse_position': f"{pair['right_primer']['start']}-{pair['right_primer']['start'] + pair['right_primer']['length']}",
                                'product_size': pair['product_size'],
                                'gc_content': (pair['left_primer']['analysis']['gc_content'] + 
                                             pair['right_primer']['analysis']['gc_content']) / 2,
                                'quality_scores': {
                                    'overall_score': pair['overall_quality'],
                                    'length_score': (pair['left_primer']['analysis']['quality_scores']['length_score'] + 
                                                   pair['right_primer']['analysis']['quality_scores']['length_score']) / 2,
                                    'gc_score': (pair['left_primer']['analysis']['quality_scores']['gc_score'] + 
                                               pair['right_primer']['analysis']['quality_scores']['gc_score']) / 2,
                                    'tm_score': (pair['left_primer']['analysis']['quality_scores']['tm_score'] + 
                                               pair['right_primer']['analysis']['quality_scores']['tm_score']) / 2,
                                    'complexity_score': (pair['left_primer']['analysis']['quality_scores']['complexity_score'] + 
                                                       pair['right_primer']['analysis']['quality_scores']['complexity_score']) / 2,
                                    'terminal_score': (pair['left_primer']['analysis']['quality_scores']['terminal_score'] + 
                                                     pair['right_primer']['analysis']['quality_scores']['terminal_score']) / 2
                                }
                            }
                            all_results.append(primer_data)
                        
                except Exception as e:
                    print(f"[PRINT-DEBUG] Error processing sequence {seq_data.get('name', 'unknown')}: {str(e)}")
                    current_app.logger.error(f"Error processing sequence {seq_data['name']}: {str(e)}")
                    continue
            
            # Save results
            if all_results:
                db_manager.save_job_results(job_id, {'primers': all_results})
                db_manager.save_primer_results(job_id, all_results)
                db_manager.update_job_status(job_id, 'completed')
            else:
                db_manager.update_job_status(job_id, 'failed', 'No primers could be designed')
                
        except Exception as e:
            print(f"[PRINT-DEBUG] Error in primer design processing for job {job_id}: {str(e)}")
            error_message = f"Error in primer design processing: {str(e)}"
            current_app.logger.error(error_message)
            db_manager.update_job_status(job_id, 'failed', error_message)

@api_bp.route('/visualize/<job_id>/<pair_index>', methods=['GET'])
def get_visualizations(job_id, pair_index):
    job = db_manager.get_job(job_id)
    if not job or not job['results']:
        return jsonify({'error': 'Job not found'}), 404
    
    try:
        pair_index = int(pair_index)
        pair = job['results']['primer_pairs'][pair_index]
        sequence = job['input_data']['sequence']  # Original sequence
        
        viz = VisualizationEngine()
        
        return jsonify({
            'binding_sites': viz.plot_binding_sites(
                sequence, 
                [pair['left_primer'], pair['right_primer']]
            ),
            'melting_profile': viz.plot_melting_profile(sequence),
            'dimer_structures': [
                viz.visualize_dimer(d) 
                for d in pair['dimer_analysis']
            ]
        })
    except Exception as e:
        return jsonify({'error': str(e)}), 500