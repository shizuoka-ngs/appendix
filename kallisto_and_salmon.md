# 発現定量化ツールsalmonによる定量とkallistoとの比較

[AJACSa6](https://github.com/AJACS-training/AJACSa6/tree/master/03_bono)の
発展課題を


## インストール

作業はconda環境で行います。必要なツールは可能な限りconda installします。

### minicondaインストール
[Condaの公式サイト](https://docs.conda.io/en/latest/miniconda.html)から
自分の環境にあったインストーラをDLしてください。

### condaのコマンド
[anaconda のコマンドリストメモ](https://qiita.com/natsuriver/items/4ae6eed5f47e34817090)

### condaで環境を作る

```
$ conda create -n eg-ngs python=3.7
$ source actibate eg-ngs
$ conda config --add channels defaults
$ conda config --add channels conda-forge
$ conda config --add channels bioconda

$ conda install -c bioconda sra-tools 
$ conda install kallisto
$ conda install salmon

$ conda install jupyter jupyter_console qtconsole notebook nbconvert
$ conda install seaborn

```

※ DRR100656,DRR100657での例。このサンプルでkallistoとsalmonの定量を行い結果を比較します。

## FASTQのトリミング

`trim_galore --fastqc --trim1 --gzip --paired DRR100657.sra_1.fastq DRR100657.sra_2.fastq`

詳細は[AJACSa6 遺伝子発現データ解析の実際](https://github.com/AJACS-training/AJACSa6/tree/master/03_bono)で

## salmonで発現定量

### salmon検索用のインデックス作成

```
$ wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.transcripts.fa.gz
$ salmon index -p 2 -t gencode.v31.transcripts.fa.gz -i gencode_v31
```

### 発現定量

```
salmon quant -p 2 -i gencode_v31 -l A --validateMappings -1 DRR100656.sra_1_val_1.fq.gz  -2 DRR100656.sra_2_val_2.fq.gz -o salmon_output_DRR100656
salmon quant -p 2 -i gencode_v31 -l A --validateMappings -1 DRR100657.sra_1_val_1.fq.gz  -2 DRR100657.sra_2_val_2.fq.gz -o salmon_output_DRR100657
```

## 可視化

### 定量結果の確認

```
$ ls salmon_output_DRR100656
aux_info/  cmd_info.json  libParams/  lib_format_counts.json  logs/  quant.sf
$ ls salmon_output_DRR100656
aux_info/  cmd_info.json  libParams/  lib_format_counts.json  logs/  quant.sf
# それぞれのディレクトリのquant.sfが発現定量結果
$ less salmon_output_DRR100656/quant.sf
$ less salmon_output_DRR100657/quant.sf
# 左から4カラム目が"TPM"値なので、高いもの順に並べて見てみる。'q'を押すと終了する
# Nameは|で繋がって妙に長い
$ sort -k 4 -rn salmon_output_DRR100656/quant.sf | less
$ sort -k 4 -rn salmon_output_DRR100657/quant.sf | less
```

以下jupyter notebookでの操作。kallistoとカラム名が微妙に違うので留意。
```
import pandas as pd
e1 = pd.read_table('salmon_output_DRR100656/quant.sf')
e1 = e1.drop(columns=['Length', 'EffectiveLength', 'NumReads'])
e1.columns = ['Name', 'TPM_DRR100656']
e2 = pd.read_table('salmon_output_DRR100657/quant.sf')
e2 = e2.drop(columns=['Length', 'EffectiveLength', 'NumReads'])
e2.columns = ['Name', 'TPM_DRR100657']

# 二つのデータを'target_id'で結合
e = pd.merge(e1, e2, on='Name')

# DataFrameに必要なフィールドが含まれていることを確認
e.head()

```

```
%matplotlib inline
import seaborn as sns

ax = sns.scatterplot(x="TPM_DRR100656", y="TPM_DRR100657", data=e)
```

jointplotだと

```
# jointplot
sns.jointplot(e[("TPM_SRR7300567")], e[("TPM_SRR7300569")])

```

TPMの対数を計算しdfに追加＆散布図にプロットする
```
import numpy as np
e['log_DRR100656'] = np.log10(e['TPM_DRR100656'] + 1)
e['log_DRR100657'] = np.log10(e['TPM_DRR100657'] + 1)
e['diff'] = abs(e['log_DRR100656'] - e['log_DRR100657'])
```

同様にseabornでプロット

```
ax = sns.scatterplot(x='log_DRR100656', y='log_DRR100657', data=e)

# 回帰直線を重ねる場合は
ax = sns.regplot(x='log_DRR100656', y='log_DRR100657', data=e, scatter_kws={'color': 'steelblue'}, line_kws={'color': 'orange'})
```

kallistoのデータでも同様に
```
ke1 = pd.read_table("DRR100656_kallisto_out/abundance.tsv")
ke1 = ke1.drop(columns=['length', 'eff_length', 'est_counts'])
ke1.columns = ['Name', 'TPM_DRR100656']
ke2 = pd.read_table("DRR100657_kallisto_out/abundance.tsv")
ke2 = ke2.drop(columns=['length', 'eff_length', 'est_counts'])
ke2.columns = ['Name', 'TPM_DRR100657']

ke = pd.merge(ke1, ke2, on='Name')
ke['log_DRR100656'] = np.log10(ke['TPM_DRR100656'] + 1)
ke['log_DRR100657'] = np.log10(ke['TPM_DRR100657'] + 1)

ax = sns.regplot(x='log_DRR100656', y='log_DRR100657', data=ke, scatter_kws={'color': 'steelblue'}, line_kws={'color': 'orange'})
```


## kallistとsalmonの比較

### 散布図を重ねてプロット

```
import matplotlib.pyplot as plt

fig = plt.figure()
ax2 = fig.add_subplot(1,1,1)
ax2.scatter(x='log_DRR100656', y='log_DRR100657', label='salmon',  data=se, color='red', marker='.', alpha=0.8)
ax2.scatter(x='log_DRR100656', y='log_DRR100657', label='kallisto', data=ke, color='steelblue', marker='x', alpha=0.8)
ax2.legend()
```

