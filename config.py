import os
from datetime import timedelta

class Config:
    # Flask settings
    SECRET_KEY = os.environ.get('SECRET_KEY') or 'dev-secret-key-change-in-production'
    
    # File upload settings
    UPLOAD_FOLDER = os.environ.get('UPLOAD_FOLDER', 'uploads')
    RESULTS_FOLDER = os.environ.get('RESULTS_FOLDER', 'results')
    MAX_CONTENT_LENGTH = 16 * 1024 * 1024  # 16MB max file size
    ALLOWED_EXTENSIONS = {'fasta', 'fa', 'fas', 'txt'}
    
    # Database settings
    DATABASE_PATH = os.environ.get('DATABASE_PATH', 'primer_design.db')
    
    # Primer design default parameters
    PRIMER_DEFAULTS = {
        'primer_min_size': 18,
        'primer_max_size': 25,
        'primer_opt_size': 20,
        'primer_min_tm': 57.0,
        'primer_max_tm': 63.0,
        'primer_opt_tm': 60.0,
        'primer_min_gc': 20.0,
        'primer_max_gc': 80.0,
        'primer_opt_gc_percent': 50.0,
        'primer_max_poly_x': 5,
        'primer_max_ns_accepted': 0,
        'product_size_range': [[75, 100], [100, 300], [300, 500], [500, 700], [700, 1000]]
    }
    
    # BLAST settings
    BLAST_DATABASE = os.environ.get('BLAST_DATABASE', 'nt')
    BLAST_EVALUE = float(os.environ.get('BLAST_EVALUE', '10.0'))
    BLAST_WORD_SIZE = int(os.environ.get('BLAST_WORD_SIZE', '7'))
    
    # Job processing
    MAX_CONCURRENT_JOBS = int(os.environ.get('MAX_CONCURRENT_JOBS', '5'))
    JOB_TIMEOUT = int(os.environ.get('JOB_TIMEOUT', '3600'))  # 1 hour
    
    # Rate limiting
    RATE_LIMIT = os.environ.get('RATE_LIMIT', "100 per hour")
    
    @staticmethod
    def init_app(app):
        # Create necessary directories
        os.makedirs(Config.UPLOAD_FOLDER, exist_ok=True)
        os.makedirs(Config.RESULTS_FOLDER, exist_ok=True)