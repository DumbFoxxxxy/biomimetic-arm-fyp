"""Project locations and Unicode-safe MuJoCo loading; no directory links."""
from .paths import WorkspacePath, PROJECT_ROOT, project_path
from .mujoco_io import load_mujoco_model, load_mujoco_xml
