@echo off
setlocal

rem One-click MuJoCo viewer for the cable-driven arm.
rem Opens the model that is already built, so it starts in about a second.
rem After changing the SolidWorks bones, run run_cable_drive.bat instead: it
rem re-exports, re-derives the wire anchors and re-optimises the guides first.

set PYTHON=D:\MujocoEnv\Scripts\python.exe
set ROOT=%~dp0
set MODEL=%ROOT%arm_cable_drive.xml

cd /d "%ROOT%"

if not exist "%PYTHON%" (
  echo MuJoCo python environment not found:
  echo   %PYTHON%
  pause
  exit /b 1
)

if not exist "%MODEL%" (
  echo Model not built yet, running the full pipeline once.
  call "%ROOT%run_cable_drive.bat"
  exit /b %errorlevel%
)

echo Opening MuJoCo viewer with %MODEL%
echo Drag the pull_ sliders in the right-hand control panel.
"%PYTHON%" "%ROOT%run_cable_drive.py"

if errorlevel 1 (
  echo.
  echo The viewer exited with an error, see the message above.
  pause
)

endlocal
