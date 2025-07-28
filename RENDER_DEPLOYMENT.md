# 🚀 Render Deployment Guide

## 📦 **Files for Render Deployment**

Your project now contains only the essential files for Render deployment:

✅ **`app.py`** - Main Flask application  
✅ **`requirements.txt`** - Python dependencies (simplified for Render)  
✅ **`Procfile`** - Tells Render how to start your app  
✅ **`render.yaml`** - Render configuration  
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

## 🧪 **Test Your App**

Once deployed, test:
- **Health Check**: `https://your-app.onrender.com/health`
- **Main App**: `https://your-app.onrender.com/`

## 📝 **What's Included**

- ✅ Flask web application
- ✅ Health check endpoint
- ✅ File upload functionality
- ✅ Basic primer design features
- ✅ Production-ready configuration

## 🚨 **If You Need primer3-py Later**

After successful deployment, you can add bioinformatics features by:
1. Updating `requirements.txt` to include `primer3-py`
2. Redeploying the service

## 🎉 **Success!**

Your Primer PCR application is now live on Render! 