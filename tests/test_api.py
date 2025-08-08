<<<<<<< HEAD
import requests
import json
import time

# Base URL for your API
BASE_URL = "http://localhost:5000/api"

def test_json_primer_design():
    """Test primer design with JSON input"""
    
    test_data = {
        "sequences": [
            {
                "name": "test_sequence_1",
                "description": "Test sequence for primer design",
                "sequence": "ATGCGATCGCATGCGATCGCATGCGATCGCATGCGATCGCATGCGATCGCATGCGATCGCATGCGATCGCATGCGATCGCATGCGATCGCATGCGATCGCATGCGATCGCATGCGATCGCATGCGATCGC",
                "length": 140
            }
        ],
        "parameters": {
            "primer_min_size": 18,
            "primer_max_size": 25,
            "primer_opt_size": 20,
            "primer_min_tm": 57.0,
            "primer_max_tm": 63.0,
            "primer_opt_tm": 60.0,
            "primer_min_gc": 40.0,
            "primer_max_gc": 60.0,
            "product_size_range": "100-300,300-500"
        }
    }
    
    print("Testing JSON primer design...")
    
    # Submit job
    response = requests.post(f"{BASE_URL}/design", json=test_data)
    print(f"Submit response: {response.status_code}")
    print(f"Response data: {response.json()}")
    
    if response.status_code == 202:
        job_id = response.json()['job_id']
        
        # Poll for results
        while True:
            status_response = requests.get(f"{BASE_URL}/status/{job_id}")
            status_data = status_response.json()
            
            print(f"Job status: {status_data['status']}")
            
            if status_data['status'] == 'completed':
                print("Results:")
                print(json.dumps(status_data['results'], indent=2))
                break
            elif status_data['status'] == 'failed':
                print(f"Job failed: {status_data.get('error_message')}")
                break
            
            time.sleep(2)

def test_file_upload():
    """Test primer design with file upload"""
    
    # Create a test FASTA file
    fasta_content = """>test_sequence_1
ATGCGATCGCATGCGATCGCATGCGATCGCATGCGATCGCATGCGATCGCATGCGATCGC
ATGCGATCGCATGCGATCGCATGCGATCGCATGCGATCGCATGCGATCGCATGCGATCGC
ATGCGATCGCATGCGATCGCATGCGATCGCATGCGATCGCATGCGATCGCATGCGATCGC

>test_sequence_2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
"""
    
    print("Testing file upload...")
    
    # Prepare files and data
    files = {'file': ('test.fasta', fasta_content, 'text/plain')}
    data = {
        'param_primer_min_size': '18',
        'param_primer_max_size': '25',
        'param_primer_opt_tm': '60.0',
        'param_primer_min_gc': '40.0',
        'param_primer_max_gc': '60.0'
    }
    
    # Submit job
    response = requests.post(f"{BASE_URL}/design", files=files, data=data)
    print(f"Submit response: {response.status_code}")
    print(f"Response data: {response.json()}")
    
    if response.status_code == 202:
        job_id = response.json()['job_id']
        
        # Poll for results
        while True:
            status_response = requests.get(f"{BASE_URL}/status/{job_id}")
            status_data = status_response.json()
            
            print(f"Job status: {status_data['status']}")
            
            if status_data['status'] == 'completed':
                print("Results found!")
                break
            elif status_data['status'] == 'failed':
                print(f"Job failed: {status_data.get('error_message')}")
                break
            
            time.sleep(2)

def test_batch_processing():
    """Test batch processing with multiple files"""
    
    # Create multiple test files
    files = []
    for i in range(3):
        fasta_content = f""">batch_sequence_{i+1}
{'ATCG' * 50}{'GCTA' * 25}{'CGAT' * 25}
"""
        files.append(('files', (f'batch_test_{i+1}.fasta', fasta_content, 'text/plain')))
    
    print("Testing batch processing...")
    
    data = {
        'param_primer_min_size': '18',
        'param_primer_max_size': '25',
        'param_primer_opt_tm': '60.0'
    }
    
    # Submit batch job
    response = requests.post(f"{BASE_URL}/batch", files=files, data=data)
    print(f"Batch submit response: {response.status_code}")
    print(f"Response data: {response.json()}")
    
    if response.status_code == 202:
        batch_id = response.json()['batch_id']
        
        # Poll batch status
        while True:
            status_response = requests.get(f"{BASE_URL}/batch/status/{batch_id}")
            status_data = status_response.json()
            
            print(f"Batch status: {status_data['overall_status']}")
            print(f"Completed: {status_data['completed']}/{status_data['total_jobs']}")
            
            if status_data['overall_status'] in ['completed', 'partial_failure']:
                print("Batch processing finished!")
                print(f"Final status: {status_data}")
                break
            
            time.sleep(3)

def test_job_history():
    """Test job history endpoint"""
    
    print("Testing job history...")
    
    response = requests.get(f"{BASE_URL}/history?limit=10")
    print(f"History response: {response.status_code}")
    
    if response.status_code == 200:
        history_data = response.json()
        print(f"Found {len(history_data['jobs'])} recent jobs:")
        for job in history_data['jobs']:
            print(f"  - {job['id']}: {job['status']} ({job['job_type']})")

def test_health_check():
    """Test API health check"""
    
    print("Testing health check...")
    
    response = requests.get("http://localhost:5000/health")
    print(f"Health check response: {response.status_code}")
    print(f"Response data: {response.json()}")

if __name__ == "__main__":
    print("Starting API tests...\n")
    
    # Run tests
    test_health_check()
    print("\n" + "="*50 + "\n")
    
    test_json_primer_design()
    print("\n" + "="*50 + "\n")
    
    test_file_upload()
    print("\n" + "="*50 + "\n")
    
    test_batch_processing()
    print("\n" + "="*50 + "\n")
    
    test_job_history()
    print("\n" + "="*50 + "\n")
    
=======
import requests
import json
import time

# Base URL for your API
BASE_URL = "http://localhost:5000/api"

def test_json_primer_design():
    """Test primer design with JSON input"""
    
    test_data = {
        "sequences": [
            {
                "name": "test_sequence_1",
                "description": "Test sequence for primer design",
                "sequence": "ATGCGATCGCATGCGATCGCATGCGATCGCATGCGATCGCATGCGATCGCATGCGATCGCATGCGATCGCATGCGATCGCATGCGATCGCATGCGATCGCATGCGATCGCATGCGATCGCATGCGATCGC",
                "length": 140
            }
        ],
        "parameters": {
            "primer_min_size": 18,
            "primer_max_size": 25,
            "primer_opt_size": 20,
            "primer_min_tm": 57.0,
            "primer_max_tm": 63.0,
            "primer_opt_tm": 60.0,
            "primer_min_gc": 40.0,
            "primer_max_gc": 60.0,
            "product_size_range": "100-300,300-500"
        }
    }
    
    print("Testing JSON primer design...")
    
    # Submit job
    response = requests.post(f"{BASE_URL}/design", json=test_data)
    print(f"Submit response: {response.status_code}")
    print(f"Response data: {response.json()}")
    
    if response.status_code == 202:
        job_id = response.json()['job_id']
        
        # Poll for results
        while True:
            status_response = requests.get(f"{BASE_URL}/status/{job_id}")
            status_data = status_response.json()
            
            print(f"Job status: {status_data['status']}")
            
            if status_data['status'] == 'completed':
                print("Results:")
                print(json.dumps(status_data['results'], indent=2))
                break
            elif status_data['status'] == 'failed':
                print(f"Job failed: {status_data.get('error_message')}")
                break
            
            time.sleep(2)

def test_file_upload():
    """Test primer design with file upload"""
    
    # Create a test FASTA file
    fasta_content = """>test_sequence_1
ATGCGATCGCATGCGATCGCATGCGATCGCATGCGATCGCATGCGATCGCATGCGATCGC
ATGCGATCGCATGCGATCGCATGCGATCGCATGCGATCGCATGCGATCGCATGCGATCGC
ATGCGATCGCATGCGATCGCATGCGATCGCATGCGATCGCATGCGATCGCATGCGATCGC

>test_sequence_2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
"""
    
    print("Testing file upload...")
    
    # Prepare files and data
    files = {'file': ('test.fasta', fasta_content, 'text/plain')}
    data = {
        'param_primer_min_size': '18',
        'param_primer_max_size': '25',
        'param_primer_opt_tm': '60.0',
        'param_primer_min_gc': '40.0',
        'param_primer_max_gc': '60.0'
    }
    
    # Submit job
    response = requests.post(f"{BASE_URL}/design", files=files, data=data)
    print(f"Submit response: {response.status_code}")
    print(f"Response data: {response.json()}")
    
    if response.status_code == 202:
        job_id = response.json()['job_id']
        
        # Poll for results
        while True:
            status_response = requests.get(f"{BASE_URL}/status/{job_id}")
            status_data = status_response.json()
            
            print(f"Job status: {status_data['status']}")
            
            if status_data['status'] == 'completed':
                print("Results found!")
                break
            elif status_data['status'] == 'failed':
                print(f"Job failed: {status_data.get('error_message')}")
                break
            
            time.sleep(2)

def test_batch_processing():
    """Test batch processing with multiple files"""
    
    # Create multiple test files
    files = []
    for i in range(3):
        fasta_content = f""">batch_sequence_{i+1}
{'ATCG' * 50}{'GCTA' * 25}{'CGAT' * 25}
"""
        files.append(('files', (f'batch_test_{i+1}.fasta', fasta_content, 'text/plain')))
    
    print("Testing batch processing...")
    
    data = {
        'param_primer_min_size': '18',
        'param_primer_max_size': '25',
        'param_primer_opt_tm': '60.0'
    }
    
    # Submit batch job
    response = requests.post(f"{BASE_URL}/batch", files=files, data=data)
    print(f"Batch submit response: {response.status_code}")
    print(f"Response data: {response.json()}")
    
    if response.status_code == 202:
        batch_id = response.json()['batch_id']
        
        # Poll batch status
        while True:
            status_response = requests.get(f"{BASE_URL}/batch/status/{batch_id}")
            status_data = status_response.json()
            
            print(f"Batch status: {status_data['overall_status']}")
            print(f"Completed: {status_data['completed']}/{status_data['total_jobs']}")
            
            if status_data['overall_status'] in ['completed', 'partial_failure']:
                print("Batch processing finished!")
                print(f"Final status: {status_data}")
                break
            
            time.sleep(3)

def test_job_history():
    """Test job history endpoint"""
    
    print("Testing job history...")
    
    response = requests.get(f"{BASE_URL}/history?limit=10")
    print(f"History response: {response.status_code}")
    
    if response.status_code == 200:
        history_data = response.json()
        print(f"Found {len(history_data['jobs'])} recent jobs:")
        for job in history_data['jobs']:
            print(f"  - {job['id']}: {job['status']} ({job['job_type']})")

def test_health_check():
    """Test API health check"""
    
    print("Testing health check...")
    
    response = requests.get("http://localhost:5000/health")
    print(f"Health check response: {response.status_code}")
    print(f"Response data: {response.json()}")

if __name__ == "__main__":
    print("Starting API tests...\n")
    
    # Run tests
    test_health_check()
    print("\n" + "="*50 + "\n")
    
    test_json_primer_design()
    print("\n" + "="*50 + "\n")
    
    test_file_upload()
    print("\n" + "="*50 + "\n")
    
    test_batch_processing()
    print("\n" + "="*50 + "\n")
    
    test_job_history()
    print("\n" + "="*50 + "\n")
    
>>>>>>> 75f1c6a80d6b193bf8a1fa04c4fb03060d16c47f
    print("API tests completed!")