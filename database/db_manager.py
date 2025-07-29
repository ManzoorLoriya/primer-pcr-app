import sqlite3
import json
from datetime import datetime
import threading
from config import Config

class DatabaseManager:
    def __init__(self):
        self.db_path = Config.DATABASE_PATH
        self.lock = threading.Lock()
    
    def init_db(self):
        """Initialize the database with required tables"""
        with self.lock:
            conn = sqlite3.connect(self.db_path)
            cursor = conn.cursor()
            
            # Jobs table
            cursor.execute('''
                CREATE TABLE IF NOT EXISTS jobs (
                    id TEXT PRIMARY KEY,
                    status TEXT NOT NULL,
                    job_type TEXT NOT NULL,
                    input_data TEXT,
                    parameters TEXT,
                    results TEXT,
                    error_message TEXT,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    completed_at TIMESTAMP
                )
            ''')
            
            # Primers table
            cursor.execute('''
                CREATE TABLE IF NOT EXISTS primers (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    job_id TEXT NOT NULL,
                    sequence_name TEXT,
                    forward_primer TEXT,
                    reverse_primer TEXT,
                    forward_tm REAL,
                    reverse_tm REAL,
                    product_size INTEGER,
                    gc_content REAL,
                    specificity_score REAL,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    FOREIGN KEY (job_id) REFERENCES jobs (id)
                )
            ''')
            
            # BLAST results table
            cursor.execute('''
                CREATE TABLE IF NOT EXISTS blast_results (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    primer_id INTEGER,
                    target_sequence TEXT,
                    alignment_length INTEGER,
                    mismatches INTEGER,
                    gaps INTEGER,
                    e_value REAL,
                    bit_score REAL,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    FOREIGN KEY (primer_id) REFERENCES primers (id)
                )
            ''')
            
            conn.commit()
            conn.close()
    
    def create_job(self, job_id, job_type, input_data, parameters):
        """Create a new job entry"""
        with self.lock:
            conn = sqlite3.connect(self.db_path)
            cursor = conn.cursor()
            
            cursor.execute('''
                INSERT INTO jobs (id, status, job_type, input_data, parameters)
                VALUES (?, ?, ?, ?, ?)
            ''', (job_id, 'pending', job_type, json.dumps(input_data), json.dumps(parameters)))
            
            conn.commit()
            conn.close()
    
    def update_job_status(self, job_id, status, error_message=None):
        """Update job status"""
        with self.lock:
            conn = sqlite3.connect(self.db_path)
            cursor = conn.cursor()
            
            update_query = '''
                UPDATE jobs 
                SET status = ?, updated_at = CURRENT_TIMESTAMP
            '''
            params = [status]
            
            if error_message:
                update_query += ', error_message = ?'
                params.append(error_message)
            
            if status == 'completed':
                update_query += ', completed_at = CURRENT_TIMESTAMP'
            
            update_query += ' WHERE id = ?'
            params.append(job_id)
            
            cursor.execute(update_query, params)
            conn.commit()
            conn.close()
    
    def save_job_results(self, job_id, results):
        """Save job results"""
        with self.lock:
            conn = sqlite3.connect(self.db_path)
            cursor = conn.cursor()
            
            cursor.execute('''
                UPDATE jobs 
                SET results = ?, status = 'completed', completed_at = CURRENT_TIMESTAMP
                WHERE id = ?
            ''', (json.dumps(results), job_id))
            
            conn.commit()
            conn.close()
    
    def get_job(self, job_id):
        """Retrieve job information"""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        cursor.execute('''
            SELECT id, status, job_type, input_data, parameters, results, 
                   error_message, created_at, updated_at, completed_at
            FROM jobs WHERE id = ?
        ''', (job_id,))
        
        row = cursor.fetchone()
        conn.close()
        
        if row:
            return {
                'id': row[0],
                'status': row[1],
                'job_type': row[2],
                'input_data': json.loads(row[3]) if row[3] else None,
                'parameters': json.loads(row[4]) if row[4] else None,
                'results': json.loads(row[5]) if row[5] else None,
                'error_message': row[6],
                'created_at': row[7],
                'updated_at': row[8],
                'completed_at': row[9]
            }
        return None
    
    def save_primer_results(self, job_id, primers_data):
        """Save primer design results"""
        with self.lock:
            conn = sqlite3.connect(self.db_path)
            cursor = conn.cursor()
            
            for primer_data in primers_data:
                cursor.execute('''
                    INSERT INTO primers 
                    (job_id, sequence_name, forward_primer, reverse_primer, 
                     forward_tm, reverse_tm, product_size, gc_content, specificity_score)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                ''', (
                    job_id,
                    primer_data.get('sequence_name'),
                    primer_data.get('forward_primer'),
                    primer_data.get('reverse_primer'),
                    primer_data.get('forward_tm'),
                    primer_data.get('reverse_tm'),
                    primer_data.get('product_size'),
                    primer_data.get('gc_content'),
                    primer_data.get('specificity_score', 0.0)
                ))
            
            conn.commit()
            conn.close()
    
    def get_job_history(self, limit=50):
        """Get recent job history"""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        cursor.execute('''
            SELECT id, status, job_type, created_at, completed_at
            FROM jobs 
            ORDER BY created_at DESC 
            LIMIT ?
        ''', (limit,))
        
        rows = cursor.fetchall()
        conn.close()
        
        return [{'id': row[0], 'status': row[1], 'job_type': row[2], 
                'created_at': row[3], 'completed_at': row[4]} for row in rows]