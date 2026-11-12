# 交互模型查看与部署

仓库：https://github.com/DumbFoxxxxy/biomimetic-arm-fyp （公开）。

网页地址：https://dumbfoxxxxy.github.io/biomimetic-arm-fyp/ （GitHub Pages 首次部署成功后生效）。

[FYP 展示页](https://dumbfoxxxxy.github.io/biomimetic-arm-fyp/showcase/) 提供英文项目概览、MuJoCo 双视角组合动画和两个 360° 模型入口。README 内嵌同一轨迹的 GIF。展示动画为记录轨迹回放，未渲染受载 TPU 变形。其原始数据及模型指纹见 [source.json](showcase/media/source.json)。

README 使用实际模型截图作入口；点击进入独立网页后，可以拖动旋转、缩放、平移，切换骨骼、肌肉、钢丝和连接件，或开启自动旋转。完整自定义网页运行在 GitHub Pages 中。

## 模型与能力边界

默认展示原版仿生肘部的零位装配，共 96 个对象。S1-Y2 肩部通过 `?model=shoulder` 单独展示，共 191 个对象。

这不是整臂运动仿真，也不计算 TPU 变形、接触、强度或实际驱动范围。桌面 MuJoCo 查看方式见 [运行说明](../RUNNING.md)。

## 本地打开

在项目根目录运行：

```powershell
python -m http.server 8765 --directory docs
```

浏览器打开 http://127.0.0.1:8765/ 。不要直接双击 HTML；模型数据需要通过本地 HTTP 服务加载。

展示页在 http://127.0.0.1:8765/showcase/ 。如需重新渲染展示素材，在项目 Python 环境安装 Pillow 和 imageio-ffmpeg 后运行 `python tools/build_fyp_showcase_media.py`。

## GitHub Pages

提交并推送后，在仓库 Settings → Pages 中选择 Deploy from a branch，分支 main，目录 /docs，然后 Save。等待 Pages 部署成功再打开上述网页地址。若已配置此来源，后续推送会自动触发部署。

GitHub Pages 不提供 Git LFS 文件。网页模型 scene.bin 和预览图 preview.png 因此以普通 Git 文件保存；制造文件继续保留原 LFS 设置。

## 复现与检查

更新模型数据：

```powershell
python tools/build_web_model.py
```

在装有项目几何依赖的 Python 环境运行导出器。它读取原版肘部与 S1-Y2 肩部各自的装配源，校验网格文件哈希，并为两部分分别生成清单和模型数据。肘部钢丝根据原版走线参数显示。

浏览器交互检查（需要 Playwright 和脚本指定的本机 Chrome）：

```powershell
python tools/check_web_viewer.py
```

已有浏览器检查与移动端截图见 [web_viewer](../reviews/S1_2026-11-03/web_viewer/browser_checks.json)。本地全部依赖已随网页保存，加载时无需外部 CDN。渲染使用 Three.js 0.180.0；许可及下载来源校验值见 viewer/vendor/LICENSE 和 viewer/vendor/manifest.json。

参考：[GitHub Pages 设置](https://docs.github.com/en/pages/getting-started-with-github-pages/creating-a-github-pages-site)、[Git LFS 与 Pages 限制](https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-git-large-file-storage)。
