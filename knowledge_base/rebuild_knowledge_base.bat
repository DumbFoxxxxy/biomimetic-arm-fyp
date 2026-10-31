@echo off
setlocal
set PYTHONIOENCODING=utf-8
set PY=D:\MujocoEnv\Scripts\python.exe
set KB=%~dp0

echo [1/3] extracting text and metadata from every PDF under D:\XJTLU\PSP ...
"%PY%" "%KB%build_knowledge_base.py" || goto :fail

echo.
echo [2/3] building abstract digest ...
"%PY%" "%KB%make_digest.py" || goto :fail

echo.
echo [3/3] generating INDEX.md, notes and documents.csv ...
"%PY%" "%KB%generate_knowledge_base.py" || goto :fail

echo.
echo done. open "%KB%INDEX.md"
pause
exit /b 0

:fail
echo.
echo FAILED. see messages above.
pause
exit /b 1
