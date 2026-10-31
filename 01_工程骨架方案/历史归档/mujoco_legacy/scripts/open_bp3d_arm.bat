@echo off
setlocal

rem One-click MuJoCo viewer for the BodyParts3D arm model.
rem Opens the model that is already built, so it starts in about a second.
rem After changing anything in this pipeline, run run_bp3d_arm.bat instead.

set PYTHON=D:\MujocoEnv\Scripts\python.exe
set ROOT=%~dp0
set MODEL=%ROOT%arm_bp3d.xml

cd /d "%ROOT%"

if not exist "%PYTHON%" (
  echo MuJoCo python environment not found:
  echo   %PYTHON%
  pause
  exit /b 1
)

if not exist "%MODEL%" (
  echo Model not built yet, running the full pipeline once.
  call "%ROOT%run_bp3d_arm.bat"
  exit /b %errorlevel%
)

echo Opening MuJoCo viewer with %MODEL%
echo Drag the pull_ sliders in the right-hand control panel.
echo Positive is a pull, negative releases the wire.
"%PYTHON%" "%ROOT%run_bp3d_arm.py"

if errorlevel 1 (
  echo.
  echo The viewer exited with an error, see the message above.
  pause
)

endlocal
