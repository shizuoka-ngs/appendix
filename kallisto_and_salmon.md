# 発現定量化ツールsalmonによる定量とkallistoとの比較

※ DRR100656,DRR100657での例。このサンプルでkallistoとsalmonの定量を行い結果を比較します。

## FASTQのトリミング

`trim_galore --fastqc --trim1 --gzip --paired DRR100657.sra_1.fastq DRR100657.sra_2.fastq`

詳細は[AJACSa6 遺伝子発現データ解析の実際](https://github.com/AJACS-training/AJACSa6/tree/master/03_bono)で


## salmonをインストール

```sql
conda insatall salmon
```

## salmon検索用のインデックス作成

```
$ wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.transcripts.fa.gz
$ salmon index -p 2 -t gencode.v31.transcripts.fa.gz -i gencode_v31
```

## 発現定量

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

以下jupyter notebookでの操作

```
import pandas as pd
e1 = pd.read_table('salmon_output_DRR100656/quqnt.sf')
e1 = e1.drop(columns=['Length', 'EffectiveLength', 'NumReads'])
e1.columns = ['Name', 'TPM_SRR7300567']
e2 = pd.read_table('SRR7300569のabundance.tsv')
e2 = e2.drop(columns=['Length', 'EffectiveLength', 'NumReads'])
e2.columns = ['Name', 'TPM_SRR7300569']

# 二つのデータを'target_id'で結合
e = pd.merge(e1, e2, on='Name')

# DataFrameに必要なフィールドが含まれていることを確認
e.head()

`