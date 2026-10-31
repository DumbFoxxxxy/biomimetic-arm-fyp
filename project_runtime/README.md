# 无目录联接的路径与场景运行层

`paths.py` 提供WorkspacePath：保留pathlib接口，将项目历史地址在构造和连接时转换为真实目录。
历史清单不必重写，已有证据哈希保持有效；新输出使用真实相对路径。支持当前CAD Python3.10和仿真Python3.14。

`mujoco_io.py` 使用Python读取UTF-8 XML、包含文件和网格资源，以ASCII资源键传给MuJoCo内存VFS。
不改网格坐标、碰撞设置、物理参数或标准库；显式调用load_mujoco_model/load_mujoco_xml。
资源文件缺失会明确报错。包含文件递归展开并检查循环引用。

源脚本通过自身位置寻找本目录，不依赖当前工作目录或PYTHONPATH配置。
Windows创建目录联接的权限不再是运行前提。完整使用步骤见项目根目录RUNNING.md。
