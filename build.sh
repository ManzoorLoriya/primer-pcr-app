#!/bin/bash

# Build script for Render deployment
# This script installs build tools before installing Python packages

set -e

echo "🔧 Installing system dependencies..."
apt-get update
apt-get install -y gcc g++ make cmake build-essential python3-dev

echo "📦 Installing Python dependencies..."
pip install --upgrade pip setuptools wheel
pip install cython numpy

echo "🔧 Installing primer3-py with proper build tools..."
pip install primer3-py==0.6.1 --no-cache-dir

echo "📦 Installing remaining dependencies..."
pip install -r requirements.txt

echo "✅ Build completed successfully!" 