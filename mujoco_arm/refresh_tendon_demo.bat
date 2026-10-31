@echo off
setlocal
cd /d "%~dp0"
powershell -ExecutionPolicy Bypass -File export_tendon_demo_meshes.ps1
if errorlevel 1 exit /b 1
D:\MujocoEnv\Scripts\python.exe prepare_tendon_demo_meshes.py
if errorlevel 1 exit /b 1
D:\MujocoEnv\Scripts\python.exe build_tendon_demo_model.py
if errorlevel 1 exit /b 1
D:\MujocoEnv\Scripts\python.exe check_tendon_demo.py
