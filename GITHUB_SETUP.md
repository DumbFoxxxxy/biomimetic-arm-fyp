# GitHub 项目备份与使用

仓库：[DumbFoxxxxy/biomimetic-arm-fyp](https://github.com/DumbFoxxxxy/biomimetic-arm-fyp)，私有仓库。

## 上传范围

包含两条设计路线的设计文档、源码、原始骨骼与已筛选肌肉、CAD与网格输出、试验准备、文献及历史记录。
模型、CAD、图片及PDF通过Git LFS保存。两条方案的文件各保留一份真实内容，兼容目录联接不重复上传。

排除临时依赖和缓存：`tmp/`、Python缓存、知识库临时数据、编辑器临时文件、环境文件，以及本机目录联接。
`tmp/` 中的布尔/CAD vendor依赖未上传；重新生成模型前需按脚本要求安装相应依赖。现有已导出模型仍保留。

## Windows 克隆与目录入口恢复

安装Git、Git LFS及所需Python环境，克隆后执行：

```powershell
git lfs install
git lfs pull
powershell -ExecutionPolicy Bypass -File tools/restore_workspace_links.ps1
```

恢复脚本只在项目内部建立目录联接，不覆盖普通目录；创建联接可能需要系统授权。
GitHub网页只展示真实目录，原路径兼容入口由脚本恢复。请先读根目录README选择设计路线。

历史脚本和生成场景包含 `D:/XJTLU/PSP`、本机Python及外部BP3D源库路径。
克隆不是完全自动部署：换机器或换位置需核对这些路径、Python依赖和外部完整BP3D数据。
当前MuJoCo查看应使用提供的脚本，避免直接以中文真实目录路径加载XML。

## 状态与数据边界

当前肘部为设计候选，实际TPU驱动范围未知，实物试验未完成。
文献和人体模型的原有来源、授权条件保留，不在此为第三方素材另行授予许可。
本项目未自动添加开源许可证；私有仓库不会将素材设为公开。
Git LFS使用账号存储和流量配额；未修改账号套餐、预算或付款设置。
