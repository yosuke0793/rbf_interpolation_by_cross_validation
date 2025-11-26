# differential_evolution

## ディレクトリ一覧
```
.
├── 00_code
├── 01_sample
│   ├── datamake
│   └── training
```

|ディレクトリ|解説|
|:---|:---|
|00_code| rbf補間の最適化プログラムのソースコード |
|01_sample| 弾性均質体を対象とした数値計算のサンプル |

## 各ディレクトリの内容

## 00_code

### 構成
```
00_code
├── Makefile
├── main00
│   └── main00.f90
├── module
│   ├── commonval.f90
│   └── fileinfo.f90
└── utilts
    ├── util_common.f90
    ├── util_de.f90
    ├── util_math.f90
    └── util_rbf.f90
```

|ファイル|概要|
|:---|:---|
|Makefile|コンパイル用Makefile|
|main00/main.f90|メインプログラム|
|module/commonval.f90|計算に用いるグローバル変数のモジュール|
|module/fileinfo.f90|データI/Oに用いるグローバル変数のモジュール|
|utilts/util_common.f90|データI/Oや汎用的な数値計算のサブルーチン|
|util_de.f90|差分進化の計算に用いるサブルーチン|
|util_math.f90|数学的処理の計算に用いるサブルーチン|
|util_rbf.f90|RBF補間の計算に用いるサブルーチン|

## 01_sample

- 弾性材料を対象とした数値計算例

```
01_sample
├── datamake
└── training
```

|ディレクトリ|解説|
|:---:|:---:|
|datamake|教師データを作成するディレクトリ|
|training|rbf補間の学習を実施するディレクトリ|

### datamake

- 教師データを作成するディレクトリ

#### ディレクトリ構成

```
datamake
├── condition_datamake.txt
├── figures
│   └── stress_strain_plots.png
├── results
│   └── dataset.txt
└── virtual_nmt.py
```

|ファイル|概要|
|:---|:---|
|./condition_datamake.txt|`virtual_nmt.py`の実行条件ファイル|
|figures/stress_strain_plots.png|生成された教師データの応力ーひずみ関係を示す散布図|
|results/dataset.txt|`virtual_nmt.py`によって生成された教師データ|
|./virtual_nmt.py|均質弾性体の数値材料試験結果を仮想的に算出するためのpythonプログラム|

### training

- rbf補間の学習を実施するディレクトリ

#### ディレクトリ構成

```
training
├── figures
│   └── visualize_results.png
├── inputs
│   ├── condition_training.txt
│   └── dataset.txt
├── results
│   ├── best.txt
│   ├── errors.txt
│   ├── log.txt
│   ├── rbf.txt
│   ├── rbf.bin
│   ├── reprod.txt
│   ├── stplog.txt
│   └── history.txt
├── go.sh
└── visualize_results.py
```

|ファイル|概要|
|:---|:---|
|./figures/visualize_results.png|構築された代理モデルと教師データの応答を比較した図|
|./inputs/condition_training.txt|RBF補間の最適化の実施条件|
|./inputs/dataset.txt|教師データ，`rbf_interpolation/01_sample/datamake/results/`を複製|
|./results/best.txt|最適化によって得られたパラメータと最小誤差|
|./results/errors.txt|RBF補間と教師データの応答の誤差|
|./results/log.txt|計算ログ出力用ファイル（デバッグ用）|
|./results/rbf.txt|最適化によって得られた代理モデルのモデル情報（ASCII形式）|
|./results/rbf.bin|`./results/rbf.txt`のバイナリバージョン．データ配置順は同一．|
|./results/reprod.txt|RBF補間と教師データの応答|
|./results/stplog.txt|差分進化のログ（毎世代出力）|
|./results/history.txt|差分進化のログ（パラメータ更新時に出力）|
|./go.sh|実行用スクリプト|
|./visualize_results.py|`./figures/visualize_results.png`を描画するプログラム|

## 各種ファイルのフォーマット

### datamake/condition_datamake.txt

<table>
  <tr>
    <td>npath = 1000</td>
  </tr>
  <tr>
    <td>nstep = 3</td>
  </tr>
  <tr>
    <td>yng = 2.0E+05</td>
  </tr>
  <tr>
    <td>poi = 2.0E+05</td>
  </tr>
  <tr>
    <td>max_strain1 = 1.0E-2</td>
  </tr>
  <tr>
    <td>max_strain2 = 1.0E-2</td>
  </tr>
  <tr>
    <td>max_strain3 = 1.0E-2</td>
  </tr>
  <tr>
    <td>max_strain4 = 1.0E-2</td>
  </tr>
  <tr>
    <td>max_strain5 = 1.0E-2</td>
  </tr>
  <tr>
    <td>max_strain6 = 1.0E-2</td>
  </tr>
  <tr>
    <td>min_strain1 = -1.0E-2</td>
  </tr>
  <tr>
    <td>min_strain2 = -1.0E-2</td>
  </tr>
  <tr>
    <td>min_strain3 = -1.0E-2</td>
  </tr>
  <tr>
    <td>min_strain4 = -1.0E-2</td>
  </tr>
  <tr>
    <td>min_strain5 = -1.0E-2</td>
  </tr>
  <tr>
    <td>min_strain6 = -1.0E-2</td>
  </tr>
</table>

### training/inputs/condition_training.txt（Differential Evolutionの入力ファイル）

<table>
  <tdead>
  </tdead>
  <tbody>
    <tr>
      <td>/CONDT/</td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <td>flag_strategy</td>
      <td>tolerance</td>
      <td></td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <td>nopt</td>
      <td>nagent_half</td>
      <td>ncv</td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <td>mf</td>
      <td>sf</td>
      <td>mcr</td>
      <td>scr</td>
      <td></td>
    </tr>
    <tr>
      <td>/NDATA/</td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <td>nexp</td>
      <td>nobj</td>
      <td></td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <td>ndata</td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <td>/PARAM/</td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <td>npara</td>
      <td>np_fix</td>
      <td></td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <td>flag_fixed</td>
      <td>flag_expnd</td>
      <td>flag_intgr</td>
      <td>maxpara</td>
      <td>minpara</td>
    </tr>
    <tr>
      <td>flag_fixed</td>
      <td>flag_expnd</td>
      <td>flag_intgr</td>
      <td>maxpara</td>
      <td>minpara</td>
    </tr>
    <tr>
      <td>flag_fixed</td>
      <td>flag_expnd</td>
      <td>flag_intgr</td>
      <td>maxpara</td>
      <td>minpara</td>
    </tr>
  </tbody>
</table>

- 変数の概要
<table>
  <tdead>
    <tr>
      <th>変数</th>
      <th>意味</th>
      <th colspan='2'>備考</th>
    </tr>
  </tdead>
  <tbody>
    <tr>
      <td rowspan="5"> flag_strategy</td>
      <td rowspan="5">DEの世代更新戦略の制御フラグ<br>1~5を入力</td>
      <td>
          1
      </td>
      <td>
          DE/RAND/1 
      </td>
    </tr>
    <tr>
      <td>
          2
      </td>
      <td>
          DE/BEST/1
      </td>
    </tr>
    <tr>
      <td>
          3
      </td>
      <td>
          DE/CURRENT-TO/1
      </td>
    </tr>
    <tr>
      <td>
          4
      </td>
      <td>
           DE/CURRENT-TO-BEST/1
      </td>
    </tr>
    <tr>
      <td>
          5
      </td>
      <td>
          【推奨】DE/RAND-TO-BEST
      </td>
    </tr>
    <tr>
      <td>tolerance</td>
      <td>DEの収束閾値<br> 収束判定：<br>パラメータ空間中の全個体の分散< tolerance </td>
      <td colspan='2'>1.0E-3 以下を推奨</td>
    </tr>
    <tr>
      <td>nopt</td>
      <td>DEの最大更新回数（世代数）</td>
      <td colspan='2'> </td>
    </tr>
    <tr>
      <td>nagent</td>
      <td>1世代中の個体数の半分</td>
      <td colspan='2'></td>
    </tr>
    <tr>
      <td>ncv</td>
      <td>k-fold Cross-validationのバッチ数（k）</td>
      <td colspan='2'>5~10を推奨．<br>データ数を割り切れる数にする必要はない</td>
    </tr>
    <tr>
      <td>mf</td>
      <td>突然変異係数の平均</td>
      <td colspan='2'>推奨：0.5</td>
    </tr>
    <tr>
      <td>sf</td>
      <td>突然変異係数の標準偏差</td>
      <td colspan='2'>推奨：0.15</td>
    </tr>
    <tr>
      <td>mcr</td>
      <td>交差率の平均</td>
      <td colspan='2'>推奨：0.25</td>
    </tr>
    <tr>
      <td>scf</td>
      <td>突然変異係数の標準偏差</td>
      <td colspan='2'>推奨：0.075</td>
    </tr>
    <tr>
      <td rowspan='2'>nexp</td>
      <td rowspan='2'>説明変数（入力変数）の数</td>
      <td>弾性材料の場合</td>
      <td>６</td>
    </tr>
    <tr>
      <td>弾塑性材料の場合</td>
      <td>12</td>
    </tr>
    <tr>
      <td>nobj</td>
      <td>目的変数（出力変数）の数</td>
      <td colspan='2'>6</td>
    </tr>
    <tr>
      <td>ndata</td>
      <td>教師データ数</td>
      <td colspan='2'>入出力のペアの数</td>
    </tr>
    <tr>
      <td>npara</td>
      <td>パラメータ数</td>
      <td colspan='2'>3</td>
    </tr>
    <tr>
      <td>np_fix</td>
      <td>固定するパラメータの数</td>
      <td colspan='2'> np_fix <= npara </td>
    </tr>
    <tr>
      <td rowspan='3'>flag_fixed</td>
      <td rowspan='3'>パラメータ固定フラグ</td>
      <td colspan='2'>該当のパラメータを固定するフラグ</td>
    </tr>
    <tr>
      <td>1</td>
      <td>固定する</td>
    </tr>
    <tr>
      <td>0</td>
      <td>固定しない（最適化の対象）</td>
    </tr>
    <tr>
      <td rowspan='3'>flag_expnd</td>
      <td rowspan='3'>探索領域自動拡大設定フラグ</td>
      <td colspan='2'>
        探査領域の端部で暫定解が発見された場合，<br>
        該当のパラメータの探査範囲を拡大するか設定するフラグ
      </td>
    </tr>
    <tr>
      <td>1</td>
      <td>拡大する</td>
    </tr>
    <tr>
      <td>0</td>
      <td>拡大しない</td>
    </tr>
    <tr>
      <td rowspan='3'>flag_intgr</td>
      <td rowspan='3'>整数化フラグ</td>
      <td colspan='2'>該当のパラメータを整数として処理するフラグ</td>      
    </tr>
    <tr>
      <td>1</td>
      <td>整数として処理</td>
    </tr>
    <tr>
      <td>0</td>
      <td>実数として処理</td>
    </tr>
    <tr>
      <td rowspan='3'>maxpara</td>
      <td rowspan='3'>探査領域の最大値</td>
      <td>1行目</td>
      <td> $n_{\psi}$, ３</td>
    </tr>
    <tr>
      <td>2行目</td>
      <td>$\log_{10}\beta$，推奨値：2 </td>
    </tr>
    <tr>
      <td>3行目</td>
      <td>$\log_{10}\eta$，推奨値：-2 </td>
    </tr>
    <tr>
      <td rowspan='3'>maxpara</td>
      <td rowspan='3'>探査領域の最小値</td>
      <td>1行目</td>
      <td> $n_{\psi}$, 1</td>
    </tr>
    <tr>
      <td>2行目</td>
      <td>$\log_{10}\beta$，推奨値：-5 </td>
    </tr>
    <tr>
      <td>3行目</td>
      <td>$\log_{10}\eta$，推奨値：-15 </td>
    </tr>
  </tbody>
</table>

### training/inputs/dataset.txt

- 教師データをASCII形式で保持
- 行数$＝$データ数（ndata）
- 列数$＝$nexp$+$nobj

| 1列目 | $\cdots$ | nexp 列目 | nexp+1列目 | $\cdots$ | nexp$+$nobj列目 |
| --- |--- |--- |--- |--- |--- |
| 1つ目の説明変数 | $\cdots$ | nexp目の説明変数 | 1つ目の目的変数 | $\cdots$ | nobj目の目的変数 |
| $\vdots$ | $\ddots$ | $\vdots$ | $\vdots$ | $\ddots$ | $\vdots$ |
| 1つ目の説明変数 | $\cdots$ | nexp目の説明変数 | 1つ目の目的変数 | $\cdots$ | nobj目の目的変数 |


### training/results/best.txt

<table>
  <tr>
    <td>Best error</td>
  </tr>
  <tr>
    <td>最適化問題の損失の最小値</td>
  </tr>
  <tr>
    <td>Parameters</td>
  </tr>
  <tr>
    <td>最適化によって得られた$n_{\psi}$</td>
  </tr>
  <tr>
    <td>最適化によって得られた$\log_{10}\beta$</td>
  </tr>
  <tr>
    <td>最適化によって得られた$\log_{10}\eta$</td>
  </tr>
</table>

### training/results/errors.txt

- 行数$=$データ数（ndata）
- 列数$=$nexp$+$nobj
- 誤差$=$代理モデルの応答$-$教師データの応答

| 1列目 | $\cdots$ | nexp 列目 | nexp+1列目 | $\cdots$ | nexp$+$nobj列目 |
| --- |--- |--- |--- |--- |--- |
| 1つ目の説明変数 | $\cdots$ | nexp目の説明変数 | 1つ目の目的変数の誤差 | $\cdots$ | nobj目の目的変数の誤差 |
| $\vdots$ | $\ddots$ | $\vdots$ | $\vdots$ | $\ddots$ | $\vdots$ |
| 1つ目の説明変数 | $\cdots$ | nexp目の説明変数 | 1つ目の目的変数の誤差 | $\cdots$ | nobj目の目的変数の誤差 |


### training/results/rbf.txt

- 最適化によって得られたRBF補間のパラメータと教師データ，重み係数を記録

<table>
  <tr>
    <td>$n_{\psi}$</td>
    <td>1（flag_normalization）</td>
    <td></td>
  </tr>
  <tr>
    <td>ndata</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>nexp</td>
    <td>nobj</td>
    <td></td>
  </tr>
  <tr>
    <td>最適化によって得られた$\beta$</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>最適化によって得られた$\eta$</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>1番目の教師データの説明変数 １</td>
    <td>...</td>
    <td>1番目の教師データの説明変数 nexp</td>
  </tr>
  <tr>
    <td>:</td>
    <td>:</td>
    <td>:</td>
  </tr>
  <tr>
    <td>ndata番目の教師データの説明変数 １</td>
    <td>...</td>
    <td>ndata番目の教師データの説明変数 nexp</td>
  </tr>
  <tr>
    <td>1番目の教師データに対応する重み係数 １</td>
    <td>...</td>
    <td>ndata番目の教師データに対応する重み係数 １</td>
  </tr>
  <tr>
    <td>:</td>
    <td>:</td>
    <td>:</td>
  </tr>
  <tr>
    <td>1番目の教師データに対応する重み係数 nobj</td>
    <td>...</td>
    <td>ndata番目の教師データに対応する重み係数 nobj</td>
  </tr>
</table>

### training/results/rbf.bin

- ` training/results/rbf.txt`と同一の配置順

### training/results/reprod.txt

- 行数$=$データ数（ndata）
- 列数$=$nexp$+$2nobj

| 1列 | $\cdots$ | nexp 列 | nexp+1列 | $\cdots$ | nexp$+$nobj列 | nexp+nobj+1列 | $\cdots$ | nexp$+$2nobj列 |
| --- |--- |--- |--- |--- |--- |--- |--- |--- |
| 第1の説明変数 | $\cdots$ | 第nexpの説明変数 | 第1の目的変数（教師データの応答） | $\cdots$ | 第nobjの目的変数（教師データの応答） | 第1の目的変数（代理モデルの応答） | $\cdots$ | 第nobjの目的変数（代理モデルの応答） |
| $\vdots$ | $\ddots$ | $\vdots$ | $\vdots$ | $\ddots$ | $\vdots$ | $\vdots$ | $\ddots$ | $\vdots$ |
| 第1の説明変数 | $\cdots$ | 第nexpの説明変数 | 第1の目的変数（教師データの応答） | $\cdots$ | 第nobjの目的変数（教師データの応答） | 第1の目的変数（代理モデルの応答） | $\cdots$ | 第nobjの目的変数（代理モデルの応答） |

### training/results/stplog.txt

|Step|Loss|Dev.Err|Dev.Para|Time[s]|H M S|
|:---|:---|:---|:---|:---|:---|
|DEの世代|その世代までの損失の最小値|その世代における誤差の標準偏差|その世代における個体の標準偏差|経過時間(秒)|経過時間(時,分,秒)|

### training/results/history.txt

|Step|ID_agent|Loss|Diff.Error|Paraemters|||
|:---|:---|:---|:---|:---|:---|:---|
|DEの世代|最も優れた個体番号|その世代までの損失の最小値|損失の最小値の変化量|暫定解$n_{\psi}$|暫定解$\log_{10}\beta$|暫定解$\log_{10}\eta$|
