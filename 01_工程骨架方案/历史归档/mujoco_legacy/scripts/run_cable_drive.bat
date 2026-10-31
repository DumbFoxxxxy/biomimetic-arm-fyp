@echo off
setlocal

set PYTHON=D:\MujocoEnv\Scripts\python.exe
set ROOT=D:\XJTLU\PSP\mujoco_arm

cd /d "%ROOT%"

if not exist "%PYTHON%" (
  echo Python environment not found: %PYTHON%
  pause
  exit /b 1
)

echo Preparing cable-drive exports...
"%PYTHON%" "%ROOT%\prepare_cable_drive_exports.py"
if errorlevel 1 goto error

echo Checking bone pose and joint frames...
"%PYTHON%" "%ROOT%\check_bone_layout.py"
if errorlevel 1 (
  echo.
  echo WARNING: the exported bones or joint frames are not usable.
  echo The viewer will still open so you can inspect the pose.
  echo.
)

echo Deriving wire anchors from the bone surfaces...
"%PYTHON%" "%ROOT%\derive_cable_anchors.py"
if errorlevel 1 goto error

echo Building cable-drive model...
"%PYTHON%" "%ROOT%\build_cable_drive_model.py"
if errorlevel 1 goto error

echo Optimizing pulley and guide points...
"%PYTHON%" "%ROOT%\optimize_guide_points.py"
if errorlevel 1 goto error

echo Building cable-drive model...
"%PYTHON%" "%ROOT%\build_cable_drive_model.py"
if errorlevel 1 goto error

echo Opening MuJoCo viewer...
"%PYTHON%" "%ROOT%\run_cable_drive.py"
goto end

:error
echo.
echo Cable-drive launch failed. Review the error above.
pause
exit /b 1

:end
endlocal
