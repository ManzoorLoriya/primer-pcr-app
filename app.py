from flask import Flask, request, jsonify, render_template
from flask_cors import CORS
from werkzeug.utils import secure_filename
import os
import uuid
import json
from datetime import datetime
import sqlite3
import threading
import logging
from api.routes import api_bp
from database.db_manager import DatabaseManager
from config import Config

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

app = Flask(__name__)

try:
    app.config.from_object(Config)
except Exception as e:
    logger.error(f"Failed to load config: {str(e)}")
    # Set default config values
    app.config.update(
        SECRET_KEY=os.urandom(24),
        UPLOAD_FOLDER='uploads',
        MAX_CONTENT_LENGTH=16 * 1024 * 1024  # 16MB max file size
    )

CORS(app)

# Create upload folder if it doesn't exist
os.makedirs(app.config.get('UPLOAD_FOLDER', 'uploads'), exist_ok=True)

# Initialize database
try:
    db_manager = DatabaseManager()
    db_manager.init_db()
    logger.info("Database initialized successfully")
except Exception as e:
    logger.error(f"Database initialization failed: {str(e)}")
    db_manager = None

# Register API blueprint
try:
    app.register_blueprint(api_bp, url_prefix='/api')
    logger.info("API blueprint registered successfully")
except Exception as e:
    logger.error(f"Failed to register API blueprint: {str(e)}")

# Global dictionary to store job statuses
job_status = {}
job_lock = threading.Lock()

@app.route('/')
def index():
    """Serve the main page"""
    try:
        return render_template('index.html')
    except Exception as e:
        logger.error(f"Error rendering index page: {str(e)}")
        return jsonify({'error': 'Internal server error'}), 500

@app.route('/health')
def health_check():
    """Health check endpoint"""
    try:
        status = {
            'status': 'healthy',
            'timestamp': datetime.now().isoformat(),
            'database': 'connected' if db_manager else 'disconnected'
        }
        return jsonify(status)
    except Exception as e:
        logger.error(f"Health check failed: {str(e)}")
        return jsonify({'status': 'unhealthy', 'error': str(e)}), 500

@app.errorhandler(404)
def not_found(error):
    return jsonify({'error': 'Not found'}), 404

@app.errorhandler(500)
def internal_error(error):
    return jsonify({'error': 'Internal server error'}), 500

if __name__ == '__main__':
    try:
        port = int(os.environ.get('PORT', 5000))
        # Use production settings for deployment
        debug_mode = os.environ.get('FLASK_ENV') == 'development'
        app.run(debug=debug_mode, host='0.0.0.0', port=port)
    except Exception as e:
        logger.error(f"Failed to start application: {str(e)}")