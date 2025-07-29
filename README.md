# Primer Design & In-Silico PCR Web Tool

A modern web application for designing PCR primers, checking their specificity, and simulating PCR reactions. Built with Flask, Biopython, and Primer3, it provides an interactive interface for molecular biologists and bioinformaticians.

---

## Features

- **Primer Design:** Automated design of optimal primer pairs using Primer3.
- **Quality Analysis:** Calculates melting temperature (Tm), GC content, and primer quality scores.
- **Specificity Checking:** BLAST-based and in-silico PCR specificity analysis against human, mouse, E. coli, or custom templates.
- **Dimer & Hairpin Detection:** Identifies potential primer-dimer and hairpin structures.
- **Visualization:** Interactive charts for binding sites, melting profiles, dimer analysis, and quality distribution.
- **Export Options:** Download results as CSV, PDF, or JSON.
- **User-Friendly UI:** Drag-and-drop sequence upload, dark mode, and step-by-step workflow.

---

## Installation

### Prerequisites

- Python 3.8+
- [pip](https://pip.pypa.io/en/stable/)
- (Optional) [Redis](https://redis.io/) if using Celery for background jobs

### Clone the Repository

```bash
git clone <your-repo-url>
cd Primer_PCR
```

### Install Dependencies

```bash
pip install -r requirements.txt
```

---

## Configuration

You can customize settings via environment variables or by editing `config.py`. Key options:

- `SECRET_KEY`: Flask secret key
- `UPLOAD_FOLDER`, `RESULTS_FOLDER`: Paths for uploads/results
- `DATABASE_PATH`: SQLite database file
- `BLAST_DATABASE`: Default BLAST database (e.g., `nt`)
- `MAX_CONCURRENT_JOBS`, `JOB_TIMEOUT`: Job processing limits

Example (optional, for development):

```bash
export FLASK_ENV=development
export SECRET_KEY=your-secret
```

---

## Running the Application

### Development Server

```bash
python app.py
```

- The app will be available at 


## Usage Guide

### 1. Input Your Sequence

- **Upload** a FASTA file or **paste** your DNA sequence directly.
- Supported formats: `.fasta`, `.fa`, `.fas`, `.txt`

### 2. Set Primer Parameters

- Adjust primer length, Tm, GC content, and product size.
- Advanced options: specificity checking (BLAST), dimer avoidance, target organism.

### 3. Run Analysis

- Click **"Design Primers"** to start.
- Progress is shown step-by-step: sequence analysis, candidate generation, specificity check, optimization.

### 4. Review & Export Results

- Examine designed primer pairs, quality scores, and visualizations.
- Export results as CSV, PDF, or JSON.

---

## Parameter Guidelines

- **Primer Length:** 18–25 bp (optimal)
- **Melting Temperature (Tm):** 55–65°C
- **GC Content:** 40–60%
- **Product Size:** 100–500 bp for efficient amplification

---

## Project Structure

```
Primer_PCR/
│
├── app.py                # Main Flask app
├── config.py             # Configuration
├── requirements.txt      # Python dependencies
├── src/                  # Core algorithms and logic
├── api/                  # API routes
├── database/             # Database manager
├── templates/            # HTML templates (UI)
├── static/               # Static files (JS, CSS)
├── data/                 # Example FASTA files
├── results/              # Output results
├── uploads/              # Uploaded files
└── tests/                # Unit tests
```

---

## Credits

Developed by **Manzoor Mohammad Loriya**  
Contact: loryamanjue786@gmail.com

---

## Troubleshooting

- Ensure all dependencies are installed (`pip install -r requirements.txt`).
- For BLAST and in-silico PCR, internet access may be required for remote database queries.


---

## Acknowledgements

- [Flask](https://flask.palletsprojects.com/)
- [Biopython](https://biopython.org/)
- [Primer3](https://primer3.org/)
- [Plotly](https://plotly.com/)
- [Bootstrap](https://getbootstrap.com/)

---
