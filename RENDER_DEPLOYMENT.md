# 🚀 Render Deployment Guide

## 📦 **Files for Render Deployment**

Your project now contains only the essential files for Render deployment:

✅ **`app.py`** - Main Flask application  
✅ **`requirements.txt`** - Python dependencies (including primer3-py)  
✅ **`build.sh`** - Build script that installs necessary tools  
✅ **`Procfile`** - Tells Render how to start your app  
✅ **`render.yaml`** - Render configuration with build script  
✅ **`config.py`** - Application configuration  
✅ **`runtime.txt`** - Python version specification  

## 🔧 **Deploy to Render**

### **Step 1: Create Render Account**
1. Go to [render.com](https://render.com)
2. Sign up with GitHub account

### **Step 2: Create Web Service**
1. Click "New +" → "Web Service"
2. Connect your GitHub repository
3. Select your repository

### **Step 3: Configure Service**
- **Name**: `primer-pcr-app` (or any name)
- **Environment**: `Python 3`
- **Branch**: `main`
- **Root Directory**: Leave empty

### **Step 4: Add Environment Variables**
In the "Environment" tab, add:
- `FLASK_ENV` = `production`
- `SECRET_KEY` = `your-secure-secret-key-here`
- `PYTHON_VERSION` = `3.11.7`

### **Step 5: Deploy**
Click "Create Web Service" and wait for deployment.

**Note**: The build script will automatically install the necessary compilation tools for `primer3-py`.

## 🧪 **Test Your App**

Once deployed, test:
- **Health Check**: `https://your-app.onrender.com/health`
- **Main App**: `https://your-app.onrender.com/`

## 📝 **What's Included**

- ✅ Flask web application
- ✅ Health check endpoint
- ✅ File upload functionality
- ✅ Full primer design features with Primer3
- ✅ Production-ready configuration

## 🚨 **If Build Still Fails**

If the build script approach doesn't work, you can:

1. **Use Railway instead** (often better for bioinformatics packages):
   - Install Railway CLI: `npm install -g @railway/cli`
   - Run: `railway login && railway init && railway up`

2. **Use Docker deployment** (more control):
   - Install Docker Desktop
   - Run: `docker build -t primer-pcr . && docker run -p 5000:5000 primer-pcr`

## 🎉 **Success!**

Your Primer PCR application is now live on Render! 