@echo off
setlocal

rem Full BodyParts3D pipeline, then the viewer.
rem
rem   collect_bp3d_arm.py            bones, joint frames, collision proxies
rem   build_bp3d_muscle_manifest.py  fused muscles and wire centre lines
rem   build_bp3d_ligament_manifest.py ligament strips and joint ranges
rem   build_bp3d_arm_model.py        arm_bp3d.xml
rem   check_bp3d_arm.py              structure, moment arms, neutral pose
rem   run_bp3d_arm.py                viewer with the pull_ sliders
rem
rem Takes about a minute because the muscle centre lines are sliced from the
rem BodyParts3D meshes every time.

set PYTHON=D:\MujocoEnv\Scripts\python.exe
set ROOT=%~dp0

cd /d "%ROOT%"

if not exist "%PYTHON%" (
  echo MuJoCo python environment not found:
  echo   %PYTHON%
  pause
  exit /b 1
)

echo [1/6] collecting BodyParts3D bones and deriving joints
"%PYTHON%" "%ROOT%collect_bp3d_arm.py" || goto :failed

echo [2/6] fusing muscles and extracting wire centre lines
"%PYTHON%" "%ROOT%build_bp3d_muscle_manifest.py" || goto :failed

echo [3/6] generating ligament strips and joint ranges
"%PYTHON%" "%ROOT%build_bp3d_ligament_manifest.py" || goto :failed

echo [4/6] building arm_bp3d.xml
"%PYTHON%" "%ROOT%build_bp3d_arm_model.py" || goto :failed

echo [5/6] checking the model
"%PYTHON%" "%ROOT%check_bp3d_arm.py" || goto :failed

echo [6/6] opening the viewer
"%PYTHON%" "%ROOT%run_bp3d_arm.py"
goto :done

:failed
echo.
echo Pipeline stopped with an error, see the message above.
pause
exit /b 1

:done
endlocal
