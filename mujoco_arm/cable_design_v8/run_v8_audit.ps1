$ErrorActionPreference = "Stop"

$python = "D:\MujocoEnv\Scripts\python.exe"

& $python build_v8_shoulder.py
& $python audit_v8_kinematics.py
& $python audit_v8_contact.py
& $python integration_v8.py
& $python export_v8_docs.py

