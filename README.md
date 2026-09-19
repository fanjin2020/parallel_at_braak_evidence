# PET A/T 与晚期 Braak 转录组分析

本仓库提供论文 *PET-defined amyloid and tau states and late-Braak transcriptomic context in Alzheimer’s disease* 的分析代码。

ADNI 部分分析 PET 定义的淀粉样蛋白/tau（A/T）状态、连续 tau 负荷与海马体积、CDR-SB 的关系，并在 A+ 人群中开展探索性纵向分析。GSE131617 部分分析晚期 Braak 病理相关的皮层差异表达和 GO 生物过程。两个数据集分别分析，没有参与者层面的匹配。

代码输出统计结果和论文图表的数据来源。投稿用的 Word 文件和最终排版图片不由本仓库生成。

## 目录

```text
scripts/                     Python 和 R 分析脚本
config/project.example.yaml  配置模板
data/contracts/              字段定义、基因符号和 GO 术语字典
requirements.txt             Python 依赖
run_analysis.sh              Linux 完整运行入口
LICENSE                      代码许可证
```

全部 Python 和 R 分析脚本放在一级 `scripts/` 目录。V3 代码更新包需要与原项目中的数据和 `data/contracts/` 配合使用。

## 数据准备

### ADNI

ADNI 数据须通过其数据申请流程获取，本仓库不提供参与者级数据。将以下文件放入 `data/private/adni/raw/`：

```text
UCBERKELEY_AMY_6MM.csv
UCBERKELEY_TAU_6MM.csv
UCSFFSX7.csv
CDR.csv
APOERES.csv
PTDEMOG.csv
UCBERKELEYFDG_8mm.csv
DXSUM.csv
ADSL.csv
DATA_DOWNLOADED_DATE.csv
```

`DATA_DOWNLOADED_DATE.csv` 用于确定日历随访机会。该文件需包含 `data_downloaded_date` 列，日期采用 `YYYY-MM-DD` 格式，并与数据下载批次一致。

### GSE131617

将表达矩阵和供者信息表放入 `data/public/gse131617/raw/`。配置模板中的文件名为：

```text
GSE131617-GPL5175_series_matrix.txt.gz
GSE131617_subject_info.xlsx
```

若实际文件名不同，在配置中填写对应路径。表达特征映射和通路分析还需要 HuEx 平台注释包，安装方式见下文。

保留 `data/contracts/` 中的字典文件。不要用其他版本的注释文件直接替换后，仍将结果视为同一次分析。

## 软件环境

既有分析记录使用 Python 3.12.13 和 statsmodels 0.14.6。当前 Python 依赖文件采用版本范围；如需精确复现，应保留实际运行环境的包版本。

在所用 Python 环境中安装依赖：

```bash
python -m pip install -r requirements.txt
```

R 核心环境如下：

| 组件 | 版本 |
|---|---|
| R | 4.4.2 |
| Bioconductor | 3.20 |
| limma | 3.62.2 |
| AnnotationDbi | 1.68.0 |
| org.Hs.eg.db | 3.20.0 |
| GO.db | 3.20.0 |
| huex10sttranscriptcluster.db | 8.8.0 |

此外需要 yaml、readxl、DBI、RSQLite 及 limma 相关依赖。使用 R 4.4.2 的 `Rscript`，在仓库根目录执行：

```bash
Rscript scripts/install_r_packages.R
Rscript scripts/check_r_environment.R
```

安装脚本仅在准备环境时使用。核心版本不匹配时，检查脚本会停止并列出差异。完整分析入口会检查包是否可加载，但不会替代这一步精确版本检查。

## 配置

以下命令均从仓库根目录运行。首次使用时复制配置模板：

```bash
cp config/project.example.yaml config/project.yaml
```

检查 `config/project.yaml` 中的数据路径和 `runtime.rscript`。相对数据路径默认以项目根目录为基准；数据位于其他位置时，也可以填写绝对路径。

`config/project.yaml` 是本地配置，不应提交到仓库。

## 运行全部分析

在已安装依赖的 Python 环境中执行：

```bash
export RSCRIPT=/实际安装位置/bin/Rscript
bash run_analysis.sh
```

将 `RSCRIPT` 替换为 R 4.4.2 的实际可执行文件路径。启动脚本默认调用 `python3`；需要更换时可设置 `PYTHON` 环境变量。

脚本先检查输入，再运行全部分析。默认使用 2,000 次参与者 bootstrap、1,000 次随访加权 bootstrap，随机种子为 `20260915`。结果写入新的时间戳目录：

```text
results/paper_run_YYYYMMDD_HHMMSS/
```

也可以直接指定运行参数：

```bash
python scripts/run_paper.py \
  --project-root . \
  --output-dir results/paper_run_v3_01 \
  --rscript /实际安装位置/bin/Rscript \
  --bootstrap 2000 \
  --selection-bootstrap 1000 \
  --seed 20260915
```

输出目录必须是新目录或空目录。重复运行时更换目录名，避免混入上一轮结果。添加 `--dry-run` 可检查输入并打印执行命令，不拟合模型；该模式不完成 R 包检查。

## 各部分分析

| 脚本 | 作用 |
|---|---|
| `prepare_adni.py` | 整理 ADNI 原始表，构建 PET 时间对齐主表 |
| `analyze_adni.py`、`adni_models.py` | 原始 A/T 横断面模型、FDG、APOE 及敏感性分析 |
| `cross_sectional.py` | A+T+ 对 A+T− 的直接比较和 A+ 内连续 tau 分析 |
| `tau_sensitivity.py`、`spline_stability.py` | 非线性、影响点和样条曲线稳定性分析 |
| `longitudinal.py` | A+ 人群 tau 与后续 CDR-SB 变化的探索性分析 |
| `clinical_tables.py` | 临床基线、纵向队列及有无随访者的描述与比较 |
| `longitudinal_sensitivity.py` | 参与者重抽样、日历机会和结局模型敏感性分析 |
| `selection_sensitivity.py` | 随访可获得性加权、权重和平衡度，以及全流程 bootstrap |
| `stage_sensitivity.py` | 临床阶段调整和 Gaussian/fractional-logit 模型比较 |
| `prepare_gse131617.R` | 整理表达矩阵、特征注释及供者信息 |
| `analyze_gse131617.R` | 供者分块差异表达及逐脑区删除分析 |
| `prepare_pathway_expression.R` | 构建 Entrez 基因级通路分析矩阵 |
| `analyze_pathways.R` | GO 排名通路分析及脑区方向一致性检查 |
| `summarize_uty_sensitivity.R` | 供者性别构成及排除 UTY 后的基因清单汇总 |
| `go_submission_policy.py` | GO 报告条目的统一筛选 |
| `compile_results_tables.py` | 整理原始模型的论文来源表，提供临床表专用入口 |
| `run_additional.py` | 按顺序运行补充临床分析 |
| `run_paper.py` | 调度完整流程 |
| `collect_paper_results.py` | 汇总图表来源文件并打包 |

UTY 脚本不重新拟合男性亚组模型。GO 报告中排除废弃条目 `GO:0090309` 时，不重新计算原检验集合的 FDR。

## 结果位置

完整运行目录内的主要结果为：

| 相对位置 | 内容 |
|---|---|
| `results/adni/` | 原始 ADNI 模型估计和样本量 |
| `results/gse131617/` | 特征级差异表达及脑区敏感性结果 |
| `results/pathways/` | 全部通路检验及筛选后的条目 |
| `results/gse_sensitivity/` | 供者性别及 UTY 汇总 |
| `results/manuscript_tables/` | 原始分析的紧凑来源表 |
| `results/extensions/` | 补充横断面、纵向、加权和临床阶段分析 |
| `paper_results/` | 论文图表相关来源表及对应清单 |
| `paper_results.zip` | 汇总结果压缩包 |

`paper_results/paper_table_sources.csv` 记录论文项目与来源文件的对应关系。当前汇总包含 59 个来源 CSV，包括完整 S2、供者年龄汇总，以及全部基线合格者的随访概率范围。

CSV 文件名中沿用的表号不一定等于最终论文表号，应按清单和分析名称选取结果。加权模型的 GEE 区间与全流程 bootstrap 区间分别保存，不能相互替代。

如分析已经完成，只需重新汇总：

```bash
python scripts/collect_paper_results.py \
  --run-dir results/paper_run_v3_01 \
  --output-dir results/paper_run_v3_01/paper_results_new
```

`--run-dir` 指向包含 `data/`、`results/` 的完整运行目录，不是解压后的分享结果包。新输出目录和同名 ZIP 不应与已有结果冲突。


