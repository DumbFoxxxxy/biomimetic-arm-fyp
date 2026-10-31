@echo off
setlocal
cd /d D:\XJTLU\PSP\mujoco_arm
D:\MujocoEnv\Scripts\python.exe build_tendon_demo_model.py
if errorlevel 1 exit /b 1
D:\MujocoEnv\Scripts\python.exe run_tendon_demo.py
