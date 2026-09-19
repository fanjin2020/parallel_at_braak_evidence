# V3 分析代码运行环境与依赖版本

整理日期：2026-09-18。

适用范围：`parallel_at_braak_evidence` 中的 16 个 Python 脚本、7 个 R 脚本及 `run_analysis.sh`。本文依据代码导入语句、`requirements.txt`、R 安装脚本和环境检查脚本整理。

## 1. 版本信息的含义

当前 Python 依赖文件指定版本范围，而不是逐包精确版本。核心 R 版本由 `check_r_environment.R` 明确指定。最新汇总结果目录没有完整的服务器环境记录，因此本文不将这些要求等同于最近一次运行的实测版本。


## 2. 语言与执行环境

| 组件 | 版本或设置 | 说明 |
|---|---|---|
| Python | 既有分析记录：3.12.13；当前代码未锁定精确版本 | 由运行脚本的 `PYTHON` 指定，默认 `python3` |
| R | 4.4.2 | R 安装与精确版本检查脚本要求 |
| Bioconductor | 3.20 | `install_r_packages.R` 指定的发行系列 |
| Bash | 未指定版本 | 执行 `run_analysis.sh`；需要支持脚本使用的 Bash 语法 |
| Conda | 未指定版本；非模型依赖 | 用于管理 Python/R 环境，不参与统计估计 |
| 操作系统 | 提供的运行记录为 Linux；发行版和内核版本未记录 | Python/R 脚本与 Bash 启动入口需区分 |
| BLAS/LAPACK | 未记录 | 数值计算底层库，建议随环境一并保存 |

运行入口默认 R 路径为 `/home/research/miniconda3/envs/paper2-r442/bin/Rscript`，可通过 `RSCRIPT` 环境变量替换。环境名称和安装路径不是版本号。

## 3. Python 第三方依赖

以下为 `requirements.txt` 中的完整直接依赖列表。

| 包 | 当前版本约束 | 用途 |
|---|---|---|
| NumPy | `>=1.26,<2.3` | 数组、矩阵运算、随机抽样、数值汇总 |
| pandas | `>=2.2,<3` | 数据整理、连接、分组和表格输出 |
| SciPy | `>=1.12,<2` | 概率分布、统计检验与数值函数 |
| statsmodels | `>=0.14,<0.15` | 回归、HC3 标准误、GEE、GLM、多重检验校正 |
| patsy | `>=0.5.6,<2` | 模型设计矩阵和样条基函数 |
| PyYAML | `>=6,<7` | 读取及写入 YAML 配置；代码导入名为 `yaml` |

上述范围是安装要求，不保证不同版本组合得到逐位相同的浮点结果。精确重现时，应保留实际运行环境的精确版本。

`argparse`、`csv`、`json`、`pathlib`、`datetime`、`collections`、`math`、`statistics`、`warnings`、`subprocess`、`shutil`、`zipfile`、`copy`、`importlib`、`sys` 等属于 Python 标准库，不单独通过 pip 安装，版本随 Python。

`adni_models`、`longitudinal`、`cross_sectional`、`tau_sensitivity`、`longitudinal_sensitivity` 和 `go_submission_policy` 等为项目内模块，不是第三方安装包。

## 4. R 依赖

### 4.1 精确指定的核心版本

| 包或运行时 | 版本 | 用途 |
|---|---|---|
| R | 4.4.2 | R 运行时 |
| limma | 3.62.2 | 供者相关结构、线性模型、经验贝叶斯差异表达及 CAMERA 排名通路分析 |
| AnnotationDbi | 1.68.0 | 注释数据库查询和标识符映射 |
| org.Hs.eg.db | 3.20.0 | 人类基因及 GO 注释 |
| GO.db | 3.20.0 | Gene Ontology 注释数据库 |
| huex10sttranscriptcluster.db | 8.8.0 | HuEx 平台 transcript-cluster 注释 |

`check_r_environment.R` 检查这些精确版本，不一致会停止。注意：`run_paper.py` 的常规预检查只检查包能否加载，不会自动执行全部精确版本比较；重现前应单独执行该检查脚本。

### 4.2 使用但未锁定精确版本的包

| 包 | 当前版本要求 | 用途 |
|---|---|---|
| yaml | 未指定精确版本 | 读取分析配置 |
| readxl | 未指定精确版本 | 读取供者信息 Excel 文件 |
| DBI | 未指定精确版本 | 数据库连接及查询接口 |
| RSQLite | 未指定精确版本 | 读取平台注释 SQLite 数据库 |
| BiocManager | 未指定精确版本 | 安装 Bioconductor 3.20 系列包；安装阶段使用 |
| statmod | 未在当前依赖检查中单列版本 | limma 的稳健经验贝叶斯步骤及相关计算所需的间接依赖，应确认已安装 |

`base`、`stats`、`utils`、`methods` 等随 R 安装，不作为单独的 CRAN 安装依赖。Bioconductor 注释包及上述包还会加载其他间接依赖；其完整版本必须从实际环境导出，不能仅凭顶层脚本推断。

安装脚本从 CRAN 安装部分包时未锁定版本，因此在不同时间重新安装，辅助包版本可能不同。

## 5. 数据与注释资源

以下资源影响结果重现，但不是 Python/R 软件包：

- `data/contracts/adni_master_columns.csv`：ADNI 数据字段定义。
- `data/contracts/gene_symbols.csv`：基因符号对应表。
- `data/contracts/go_bp_term_dictionary.csv`：GO 术语对应表。
- `huex10sttranscriptcluster.sqlite`：平台注释数据库文件。
- ADNI 原始数据下载批次、GSE131617 输入数据及分析配置。

应保留与此次分析一致的文件，不要仅更新软件包而替换注释资源。当前代码包没有 matplotlib、seaborn、python-docx 或 PDF 工具的直接依赖；独立绘图和文档制作脚本不在本清单范围内。

## 6. 在实际运行服务器导出精确版本

在运行分析的同一 Python 环境中执行以下命令。不要改用另一台电脑的环境来代表本次运行。

```bash
conda activate adbio
python --version
python -m pip freeze --all > python_packages_actual.txt
python -c "import sys, platform; print(sys.executable); print(sys.version); print(platform.platform())" > python_runtime_actual.txt
python -c "import numpy; numpy.show_config()" > numpy_build_actual.txt
conda list --explicit > conda_python_actual.txt
```

逐项打印六个直接 Python 依赖：

```bash
python -c "from importlib.metadata import version; names=['numpy','pandas','scipy','statsmodels','patsy','PyYAML']; print('\n'.join(n+'=='+version(n) for n in names))"
```

在项目根目录，用本次分析实际使用的 Rscript 执行：

```bash
RSCRIPT=/home/research/miniconda3/envs/paper2-r442/bin/Rscript
"$RSCRIPT" scripts/check_r_environment.R > r_core_versions_actual.txt
"$RSCRIPT" -e 'write.csv(as.data.frame(installed.packages()[,c("Package","Version","Built","LibPath")]), "r_packages_actual.csv", row.names=FALSE)'
"$RSCRIPT" -e 'pkgs <- c("limma","AnnotationDbi","org.Hs.eg.db","GO.db","huex10sttranscriptcluster.db","yaml","readxl","DBI","RSQLite","statmod"); invisible(lapply(pkgs, loadNamespace)); sessionInfo()' > r_session_actual.txt
bash --version > bash_version_actual.txt
uname -a > system_actual.txt
```

`r_packages_actual.csv` 包含该 R 环境中所有已安装包（含未被本项目使用的包），`r_session_actual.txt` 补充运行平台、已加载命名空间及数值库信息。

## 7. 本文依据

- `requirements.txt`
- `scripts/check_r_environment.R`
- `scripts/install_r_packages.R`
- `scripts/run_paper.py`
- `scripts/*.py`、`scripts/*.R` 中的依赖调用
- `run_analysis.sh`


