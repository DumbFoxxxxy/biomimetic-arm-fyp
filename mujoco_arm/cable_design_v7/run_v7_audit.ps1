$ErrorActionPreference = "Stop"
$root = Split-Path -Parent $MyInvocation.MyCommand.Path
$python = "D:\MujocoEnv\Scripts\python.exe"

& $python (Join-Path $root "audit_v7.py")
& $python (Join-Path $root "export_v7_docs.py")
