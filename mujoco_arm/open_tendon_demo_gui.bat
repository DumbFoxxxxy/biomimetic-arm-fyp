@echo off
setlocal
cd /d "%~dp0"
D:\MujocoEnv\Scripts\python.exe build_tendon_demo_model.py
if errorlevel 1 exit /b 1
D:\MujocoEnv\Scripts\python.exe run_tendon_demo_gui.py
