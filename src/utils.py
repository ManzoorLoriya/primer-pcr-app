import re
import os
from datetime import datetime
from Bio import SeqIO
from io import StringIO
from config import Config

def format_primer_results(primers_data):
    """Format primer results for API response"""
    formatted_results = []
    
    for primer in primers_data:
        formatted_primer = {
            'sequence_name': primer.get('sequence_name', ''),
            'forward_primer': {
                'sequence': primer.get('forward_primer', ''),
                'tm': round(primer.get('forward_tm', 0), 2),
                'gc_content': round(primer.get('forward_gc', 0), 1),
                'length': len(primer.get('forward_primer', '')),
                'start_position': primer.get('forward_start', 0),
                'end_position': primer.get('forward_end', 0)
            },
            'reverse_primer': {
                'sequence': primer.get('reverse_primer', ''),
                'tm': round(primer.get('reverse_tm', 0), 2),
                'gc_content': round(primer.get('reverse_gc', 0), 1),
                'length': len(primer.get('reverse_primer', '')),
                'start_position': primer.get('reverse_start', 0),
                'end_position': primer.get('reverse_end', 0)
            },
            'product': {
                'size': primer.get('product_size', 0),
                'sequence': primer.get('product_sequence', ''),
                'start_position': primer.get('product_start', 0),
                'end_position': primer.get('product_end', 0)
            },
            'quality_scores': {
                'overall_score': round(primer.get('overall_score', 0), 3),
                'specificity_score': round(primer.get('specificity_score', 0), 3),
                'dimer_score': round(primer.get('dimer_score', 0), 3),
                'tm_difference': abs(primer.get('forward_tm', 0) - primer.get('reverse_tm', 0))
            },
            'warnings': primer.get('warnings', []),
            'blast_results': primer.get('blast_results', [])
        }
        
        formatted_results.append(formatted_primer)
    
    return formatted_results

def calculate_gc_content(sequence):
    """Calculate GC content percentage"""
    if not sequence:
        return 0.0
    
    gc_count = sequence.upper().count('G') + sequence.upper().count('C')
    return (gc_count / len(sequence)) * 100

def reverse_complement(sequence):
    """Generate reverse complement of DNA sequence"""
    complement_map = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
        'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
        'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H',
        'H': 'D', 'V': 'B', 'N': 'N'
    }
    
    complement = ''.join([complement_map.get(base, base) for base in sequence.upper()])
    return complement[::-1]

def estimate_melting_temperature(sequence):
    """Estimate melting temperature using nearest neighbor method (simplified)"""
    if len(sequence) < 14:
        # Wallace rule for short sequences
        return (sequence.count('A') + sequence.count('T')) * 2 + \
               (sequence.count('G') + sequence.count('C')) * 4
    else:
        # Simplified nearest neighbor calculation
        # This is a basic approximation - use primer3 for accurate calculations
        gc_content = calculate_gc_content(sequence)
        length = len(sequence)
        
        # Basic formula: Tm = 81.5°C + 16.6*(log10[Na+]) + 0.41*(%GC) - 675/length
        # Assuming [Na+] = 0.05 M
        tm = 81.5 + 16.6 * (-1.3) + 0.41 * gc_content - 675 / length
        return max(tm, 0)

def check_primer_dimers(primer1, primer2, min_score=6):
    """Check for potential primer dimers (simplified)"""
    # This is a basic implementation - use more sophisticated tools for production
    dimers = []
    
    # Check for 3' complementarity
    primer1_3prime = primer1[-5:]  # Last 5 bases
    primer2_3prime_rc = reverse_complement(primer2[-5:])
    
    # Simple alignment check
    matches = sum(1 for a, b in zip(primer1_3prime, primer2_3prime_rc) if a == b)
    
    if matches >= 3:
        dimers.append({
            'type': '3_prime_complementarity',
            'score': matches,
            'risk_level': 'HIGH' if matches >= 4 else 'MEDIUM'
        })
    
    return dimers

def validate_primer_pair(forward_primer, reverse_primer, parameters):
    """Validate a primer pair against design parameters"""
    warnings = []
    
    # Check primer lengths
    if len(forward_primer) < parameters.get('min_length', 18):
        warnings.append(f"Forward primer too short ({len(forward_primer)} bp)")
    if len(reverse_primer) < parameters.get('min_length', 18):
        warnings.append(f"Reverse primer too short ({len(reverse_primer)} bp)")
    
    if len(forward_primer) > parameters.get('max_length', 25):
        warnings.append(f"Forward primer too long ({len(forward_primer)} bp)")
    if len(reverse_primer) > parameters.get('max_length', 25):
        warnings.append(f"Reverse primer too long ({len(reverse_primer)} bp)")
    
    # Check GC content
    forward_gc = calculate_gc_content(forward_primer)
    reverse_gc = calculate_gc_content(reverse_primer)
    
    min_gc = parameters.get('min_gc', 20)
    max_gc = parameters.get('max_gc', 80)
    
    if forward_gc < min_gc or forward_gc > max_gc:
        warnings.append(f"Forward primer GC content ({forward_gc:.1f}%) outside range")
    if reverse_gc < min_gc or reverse_gc > max_gc:
        warnings.append(f"Reverse primer GC content ({reverse_gc:.1f}%) outside range")
    
    # Check melting temperature
    forward_tm = estimate_melting_temperature(forward_primer)
    reverse_tm = estimate_melting_temperature(reverse_primer)
    
    min_tm = parameters.get('min_tm', 55)
    max_tm = parameters.get('max_tm', 65)
    
    if forward_tm < min_tm or forward_tm > max_tm:
        warnings.append(f"Forward primer Tm ({forward_tm:.1f}°C) outside range")
    if reverse_tm < min_tm or reverse_tm > max_tm:
        warnings.append(f"Reverse primer Tm ({reverse_tm:.1f}°C) outside range")
    
    # Check Tm difference
    tm_diff = abs(forward_tm - reverse_tm)
    max_tm_diff = parameters.get('max_tm_diff', 5)
    
    if tm_diff > max_tm_diff:
        warnings.append(f"Tm difference ({tm_diff:.1f}°C) exceeds limit")
    
    # Check for dimers
    dimers = check_primer_dimers(forward_primer, reverse_primer)
    if dimers:
        warnings.append(f"Potential dimer formation detected")
    
    return warnings

def generate_error_response(error_type, message, details=None):
    """Generate standardized error response"""
    response = {
        'success': False,
        'error_type': error_type,
        'message': message,
        'timestamp': datetime.now().isoformat()
    }
    
    if details:
        response['details'] = details
    
    return response

def sanitize_filename(filename):
    """Sanitize filename for safe storage"""
    # Remove or replace unsafe characters
    filename = re.sub(r'[<>:"/\\|?*]', '_', filename)
    return filename[:255]  # Limit length

def paginate_results(results, page=1, per_page=20):
    """Paginate results for API response"""
    total = len(results)
    start_idx = (page - 1) * per_page
    end_idx = start_idx + per_page
    
    return {
        'results': results[start_idx:end_idx],
        'pagination': {
            'page': page,
            'per_page': per_page,
            'total': total,
            'pages': (total + per_page - 1) // per_page
        }
    }

def allowed_file(filename):
    """Check if file extension is allowed"""
    allowed_extensions = {'fasta', 'fa', 'fas', 'txt'}
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in allowed_extensions

def parse_fasta_content(content):
    """Parse FASTA content from string"""
    sequences = []
    
    try:
        # Use BioPython to parse FASTA
        fasta_io = StringIO(content)
        for record in SeqIO.parse(fasta_io, "fasta"):
            sequences.append({
                'id': record.id,
                'name': record.name,
                'description': record.description,
                'sequence': str(record.seq),
                'length': len(record.seq)
            })
    except Exception as e:
        # Fallback to manual parsing
        sequences = parse_fasta_manual(content)
    
    return sequences

def parse_fasta_manual(content):
    """Manual FASTA parsing as fallback"""
    sequences = []
    current_sequence = None
    
    for line in content.split('\n'):
        line = line.strip()
        if not line:
            continue
            
        if line.startswith('>'):
            # Save previous sequence
            if current_sequence:
                sequences.append(current_sequence)
            
            # Start new sequence
            header = line[1:]  # Remove '>'
            parts = header.split(' ', 1)
            
            current_sequence = {
                'id': parts[0],
                'name': parts[0],
                'description': parts[1] if len(parts) > 1 else '',
                'sequence': '',
                'length': 0
            }
        else:
            # Add to current sequence
            if current_sequence:
                current_sequence['sequence'] += line.upper()
                current_sequence['length'] = len(current_sequence['sequence'])
    
    # Add last sequence
    if current_sequence:
        sequences.append(current_sequence)
    
    return sequences

def validate_parameters(parameters):
    """Validate primer design parameters"""
    errors = []
    
    # Check required parameters
    required_params = ['min_length', 'max_length', 'min_tm', 'max_tm']
    for param in required_params:
        if param not in parameters:
            errors.append(f"Missing required parameter: {param}")
    
    # Check parameter ranges
    if 'min_length' in parameters and 'max_length' in parameters:
        if parameters['min_length'] > parameters['max_length']:
            errors.append("min_length cannot be greater than max_length")
    
    if 'min_tm' in parameters and 'max_tm' in parameters:
        if parameters['min_tm'] > parameters['max_tm']:
            errors.append("min_tm cannot be greater than max_tm")
    
    if 'min_gc' in parameters and 'max_gc' in parameters:
        if parameters['min_gc'] > parameters['max_gc']:
            errors.append("min_gc cannot be greater than max_gc")
    
    # Check reasonable ranges
    if 'min_length' in parameters and parameters['min_length'] < 10:
        errors.append("min_length should be at least 10 bp")
    
    if 'max_length' in parameters and parameters['max_length'] > 50:
        errors.append("max_length should not exceed 50 bp")
    
    return errors

def validate_sequence(sequence):
    """Validate DNA sequence"""
    if not sequence:
        return False, "Empty sequence"
    
    # Check for valid DNA characters
    valid_bases = set('ATCGN')
    sequence_upper = sequence.upper()
    
    invalid_bases = set(sequence_upper) - valid_bases
    if invalid_bases:
        return False, f"Invalid bases found: {', '.join(invalid_bases)}"
    
    # Check minimum length
    if len(sequence) < 50:
        return False, "Sequence too short (minimum 50 bp)"
    
    # Check maximum length
    if len(sequence) > 100000:
        return False, "Sequence too long (maximum 100,000 bp)"
    
    return True, "Valid sequence"
