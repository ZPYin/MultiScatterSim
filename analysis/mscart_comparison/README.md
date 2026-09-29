# MSCART 同输入模型比对

本资料包用于 `GitHub-repository/MultiScatterSim-dev` 的正式 CPU/GPU 模型验证。



两个案例目录：

| 目录 | 保留内容 |
| --- | --- |
| [ground](ground/) | 地基案例：低、高消光各一套输入和结果，包含0.5、2.0 mrad全视场 |
| [satellite](satellite/) | Hu2001星载案例：一套输入和结果，包含0.001、0.13、13 mrad全视场 |

先看 [地基对比图](ground/comparison.png)、[星载对比图](satellite/comparison.png)，以及各目录的 `summary.csv`。PDF用于打印；`comparison.csv` 是图对应的平均双通道、逐箱和累计退偏比数据。


## 计算定义

线偏振通道：`I∥=(I+Q)/2`、`I⊥=(I−Q)/2`，`δ=I⊥/I∥`。先平均通道，再累加取比；阴影为独立重复范围。

- 地基：355 nm，5 m分箱；主表统计云内20–120 m原始信号。MSCART每km²输出除以10⁶换为每m²。图显示累计到150 m。
- 星载：532 nm，3 m分箱；主表统计云内0–150、0–300 m距离校正信号。图显示累计到300 m。
- 每次本模型运行200万光子；地基MSCART每次200万，星载每次1000万。4+阶为总信号减前三阶。

## 重绘

在 PyCharm 中打开本目录的 `compare.py`，右键选择“运行”，或执行 `python compare.py`。默认绘制地基、星载两个案例，读取各自的 `model_results.csv`。控制台会打印读取的 CSV 和图片保存位置；在 Windows 上完成后自动打开 PNG 图片。

使用 `--case ground` 或 `--case satellite` 可只画一个案例。添加 `--no-show` 可只保存文件。

```bash
python -m pip install -r requirements.txt
python compare.py --case ground
python compare.py --case satellite
```

只使用包内数据。

## 重跑本模型

在MATLAB中进入 `ground` 或 `satellite` 目录，打开该目录的 `run_model.m` 后直接点击“运行”，或在命令窗口执行：

```matlab
run_model
```

入口会根据自身文件位置自动定位本项目根目录。默认每次200万光子、3个随机种子；地基计算低/高消光各两种视场，星载计算三种视场。

也可手动指定参数；第一个参数为 `[]` 时仍自动定位项目根目录：

```matlab
run_model([],2000000,[20260929,20260930,20261001]);
```

显式传入 `run_model(modelRoot,nPhoton,seeds)` 的原用法继续支持。

每次运行直接覆盖该目录的 `model_results.csv`，只保留当前结果。从资料包根目录重绘：

```bash
python compare.py --case ground
python compare.py --case satellite
```

## 重跑MSCART

```bash
bash run_mscart.sh /absolute/path/to/MSCART ground_low /absolute/path/to/new_low_run
bash run_mscart.sh /absolute/path/to/MSCART ground_high /absolute/path/to/new_high_run
bash run_mscart.sh /absolute/path/to/MSCART satellite /absolute/path/to/new_satellite_run
```

`source_info.json` 记录数据来源、SHA256和整理核对结果。


