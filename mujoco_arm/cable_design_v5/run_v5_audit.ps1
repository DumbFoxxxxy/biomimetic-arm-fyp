$ErrorActionPreference = "Stop"
$Python = "D:\MujocoEnv\Scripts\python.exe"
$Root = Split-Path -Parent $MyInvocation.MyCommand.Path

& $Python (Join-Path $Root "build_v5.py")
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

& $Python (Join-Path $Root "audit_v5.py")
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

& $Python (Join-Path $Root "export_v5_docs.py")
exit $LASTEXITCODE
