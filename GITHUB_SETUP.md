# GitHub 备份与克隆

公开仓库：[DumbFoxxxxy/biomimetic-arm-fyp](https://github.com/DumbFoxxxxy/biomimetic-arm-fyp)。
包含两套方案、源码、设计与历史记录、原始输入、CAD和制造网格；大文件通过Git LFS保存。
排除临时依赖、缓存和环境文件。文献及模型保留原授权，本仓库未为第三方素材另授许可。

克隆后执行git lfs install与git lfs pull，再配置所需Python依赖即可使用查看脚本。
无需恢复目录联接。旧restore_workspace_links.ps1已退役，不再创建目录。
路径、依赖与验证范围详见 [RUNNING.md](RUNNING.md)。
外部完整BP3D库、CAD软件和tmp内的本机依赖未打包；导出模型已保存。
