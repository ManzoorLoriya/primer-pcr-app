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
            'type': '3prime_dimer',
            'score': matches,
            'description': f"Potential 3' dimer with {matches} matches"
        })
    
    # Check for hairpin formation (self-complementarity)
    for primer in [primer1, primer2]:
        rc_primer = reverse_complement(primer)
        max_matches = 0
        
        # Sliding window comparison
        for i in range(len(primer) - 4):
            for j in range(len(rc_primer) - 4):
                window1 = primer[i:i+5]
                window2 = rc_primer[j:j+5]
                matches = sum(1 for a, b in zip(window1, window2) if a == b)
                max_matches = max(max_matches, matches)
        
        if max_matches >= 4:
            dimers.append({
                'type': 'hairpin',
                'score': max_matches,
                'primer': primer,
                'description': f"Potential hairpin with {max_matches} matches"
            })
    
    return dimers

def validate_primer_pair(forward_primer, reverse_primer, parameters):
    """Validate a primer pair against design criteria"""
    issues = []
    warnings = []
    
    # Check Tm difference
    forward_tm = estimate_melting_temperature(forward_primer)
    reverse_tm = estimate_melting_temperature(reverse_primer)
    tm_diff = abs(forward_tm - reverse_tm)
    
    if tm_diff > 5:
        issues.append(f"Large Tm difference: {tm_diff:.1f}°C")
    elif tm_diff > 3:
        warnings.append(f"Moderate Tm difference: {tm_diff:.1f}°C")
    
    # Check GC content
    forward_gc = calculate_gc_content(forward_primer)
    reverse_gc = calculate_gc_content(reverse_primer)
    
    min_gc = parameters.get('primer_min_gc', 20)
    max_gc = parameters.get('primer_max_gc', 80)
    
    if forward_gc < min_gc or forward_gc > max_gc:
        issues.append(f"Forward primer GC content out of range: {forward_gc:.1f}%")
    
    if reverse_gc < min_gc or reverse_gc > max_gc:
        issues.append(f"Reverse primer GC content out of range: {reverse_gc:.1f}%")
    
    # Check for primer dimers
    dimers = check_primer_dimers(forward_primer, reverse_primer)
    for dimer in dimers:
        if dimer['score'] >= 4:
            issues.append(f"Potential {dimer['type']}: {dimer['description']}")
        else:
            warnings.append(f"Weak {dimer['type']}: {dimer['description']}")
    
    # Check for poly-runs
    max_poly = parameters.get('primer_max_poly_x', 5)
    
    for primer, name in [(forward_primer, 'Forward'), (reverse_primer, 'Reverse')]:
        for base in 'ATCG':
            poly_run = base * (max_poly + 1)
            if poly_run in primer:
                issues.append(f"{name} primer contains poly-{base} run")
    
    return {
        'issues': issues,
        'warnings': warnings,
        'tm_difference': tm_diff,
        'forward_tm': forward_tm,
        'reverse_tm': reverse_tm,
        'forward_gc': forward_gc,
        'reverse_gc': reverse_gc
    }

def generate_error_response(error_type, message, details=None):
    """Generate standardized error response"""
    response = {
        'error': True,
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
    filename = re.sub(r'[^\w\-_\.]', '_', filename)
    # Limit length
    if len(filename) > 100:
        name, ext = os.path.splitext(filename)
        filename = name[:95] + ext
    
    return filename

def paginate_results(results, page=1, per_page=20):
    """Paginate results list"""
    total = len(results)
    start = (page - 1) * per_page
    end = start + per_page
    
    return {
        'results': results[start:end],
        'pagination': {
            'page': page,
            'per_page': per_page,
            'total': total,
            'pages': (total + per_page - 1) // per_page,
            'has_next': end < total,
            'has_prev': page > 1
        }
    }
def allowed_file(filename):
    """Check if file extension is allowed"""
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in Config.ALLOWED_EXTENSIONS

def parse_fasta_content(content):
    """Parse FASTA content and return list of sequence dictionaries"""
    sequences = []
    
    try:
        # Handle different newline formats
        content = content.replace('\r\n', '\n').replace('\r', '\n')
        
        # Parse using Biopython
        fasta_io = StringIO(content)
        
        for record in SeqIO.parse(fasta_io, "fasta"):
            if len(record.seq) > 50:  # Minimum sequence length
                sequences.append({
                    'name': record.id,
                    'description': record.description,
                    'sequence': str(record.seq).upper(),
                    'length': len(record.seq)
                })
        
        # If Biopython parsing fails, try manual parsing
        if not sequences:
            sequences = parse_fasta_manual(content)
            
    except Exception as e:
        # Fallback to manual parsing
        sequences = parse_fasta_manual(content)
    
    return sequences

def parse_fasta_manual(content):
    """Manual FASTA parsing as fallback"""
    sequences = []
    current_header = None
    current_sequence = []
    
    for line in content.split('\n'):
        line = line.strip()
        if not line:
            continue
            
        if line.startswith('>'):
            # Save previous sequence if exists
            if current_header and current_sequence:
                seq_str = ''.join(current_sequence).upper()
                if len(seq_str) > 50:
                    sequences.append({
                        'name': current_header.split()[0],
                        'description': current_header,
                        'sequence': seq_str,
                        'length': len(seq_str)
                    })
            
            # Start new sequence
            current_header = line[1:]  # Remove '>'
            current_sequence = []
        else:
            # Add to current sequence
            # Remove any non-nucleotide characters
            clean_seq = re.sub(r'[^ATCGRYSWKMBDHVN]', '', line.upper())
            if clean_seq:
                current_sequence.append(clean_seq)
    
    # Save last sequence
    if current_header and current_sequence:
        seq_str = ''.join(current_sequence).upper()
        if len(seq_str) > 50:
            sequences.append({
                'name': current_header.split()[0],
                'description': current_header,
                'sequence': seq_str,
                'length': len(seq_str)
            })
    
    return sequences

def validate_parameters(parameters):
    """Validate and set default parameters"""
    defaults = Config.PRIMER_DEFAULTS.copy()
    
    # Parameter mapping from frontend to backend
    param_mapping = {
        'minLength': 'primer_min_size',
        'maxLength': 'primer_max_size',
        'minTm': 'primer_min_tm',
        'maxTm': 'primer_max_tm',
        'minGC': 'primer_min_gc',
        'maxGC': 'primer_max_gc',
        'minProduct': 'product_size_range',
        'maxProduct': 'product_size_range'
    }
    
    # Override with provided parameters
    for key, value in parameters.items():
        # Map frontend parameter to backend parameter
        backend_key = param_mapping.get(key, key)
        
        if backend_key in defaults:
            try:
                # Convert to appropriate type
                if isinstance(defaults[backend_key], float):
                    defaults[backend_key] = float(value)
                elif isinstance(defaults[backend_key], int):
                    defaults[backend_key] = int(value)
                elif isinstance(defaults[backend_key], list):
                    # Handle special cases like product_size_range
                    if backend_key == 'product_size_range':
                        if key == 'minProduct':
                            defaults[backend_key] = [[value, defaults[backend_key][0][1]]]
                        elif key == 'maxProduct':
                            defaults[backend_key] = [[defaults[backend_key][0][0], value]]
                        else:
                            defaults[backend_key] = value
                    else:
                        defaults[backend_key] = value
                else:
                    defaults[backend_key] = value
            except (ValueError, TypeError):
                # Keep default value if conversion fails
                pass
    
    # Validate parameter ranges
    if defaults['primer_min_size'] > defaults['primer_max_size']:
        defaults['primer_min_size'] = defaults['primer_max_size'] - 2
    
    if defaults['primer_min_tm'] > defaults['primer_max_tm']:
        defaults['primer_min_tm'] = defaults['primer_max_tm'] - 3
    
    if defaults['primer_min_gc'] > defaults['primer_max_gc']:
        defaults['primer_min_gc'] = defaults['primer_max_gc'] - 10
    
    return defaults

def validate_sequence(sequence):
    """Validate DNA sequence"""
    if not sequence:
        return False, "Empty sequence"
    
    # Check length
    if len(sequence) < 100:
        return False, "Sequence too short (minimum 100 bp)"
    
    if len(sequence) > 50000:
        return False, "Sequence too long (maximum 50,000 bp)"
    
    # Check for valid nucleotides
    valid_nucleotides = set('ATCGRYSWKMBDHVN')
    sequence_set = set(sequence.upper())
    
    if not sequence_set.issubset(valid_nucleotides):
        invalid_chars = sequence_set - valid_nucleotides
        return False, f"Invalid nucleotides found: {', '.join(invalid_chars)}"
    
    # Check for reasonable nucleotide composition
    if sequence.count('N') / len(sequence) > 0.1:
        return False, "Too many ambiguous nucleotides (N)"
    
    return True, "Valid sequence"

