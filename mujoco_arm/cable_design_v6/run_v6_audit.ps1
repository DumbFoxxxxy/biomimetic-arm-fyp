$ErrorActionPreference = "Stop"
$root = Split-Path -Parent $MyInvocation.MyCommand.Path
$python = "D:\MujocoEnv\Scripts\python.exe"

& $python (Join-Path $root "build_v6.py")
& $python (Join-Path $root "audit_v6.py")
& $python (Join-Path $root "export_v6_docs.py")
