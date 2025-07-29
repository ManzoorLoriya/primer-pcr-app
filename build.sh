#!/bin/bash

# Build script for Render deployment
# This script installs build tools before installing Python packages

set -e

echo "🔧 Installing system dependencies..."
apt-get update
apt-get install -y gcc g++ make cmake build-essential

echo "📦 Installing Python dependencies..."
pip install --upgrade pip
pip install -r requirements.txt

echo "✅ Build completed successfully!" 