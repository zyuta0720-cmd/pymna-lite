# pymna-lite
Python-based Circuit Analysis Tool (MNA Engine)

## 概要 / Overview
`pymna-lite` は、Pythonで実装された修正節点解析（MNA: Modified Nodal Analysis）エンジンを搭載した回路解析ツールです。直感的なGUIを備え、日本語と英語の両方のインターフェースに対応しています。

`pymna-lite` is a circuit analysis tool powered by a Python-implemented Modified Nodal Analysis (MNA) engine. It features an intuitive GUI and supports both Japanese and English interfaces.

## 主な機能 / Key Features
* **MNA Engine**: 抵抗(R)、独立電源(V/I)に加え、電圧制御電圧源(E)や電流制御電流源(F)の解析に対応。
* **Unit Support**: `k`, `meg`, `m`, `u`, `n`, `p` などの単位を自動で数値変換します。
* **LTspice Integration**: 解析結果を `.asc` および `.cir` ファイルとして書き出し可能です。

* **MNA Engine**: Supports analysis of Resistors (R), Independent Sources (V/I), Voltage-Controlled Voltage Sources (E), and Current-Controlled Current Sources (F).
* **Unit Support**: Automatically converts units such as `k`, `meg`, `m`, `u`, `n`, and `p` to numerical values.
* **LTspice Integration**: Export your analysis to `.asc` and `.cir` files for use in LTspice.

## 使い方 / How to Use
1. **Requirements**: `pip install numpy` で必要なライブラリをインストールしてください。
   (Install the required library via `pip install numpy`)
2. **Run**: `pymna-lite.py` を実行します。
   (Run `pymna-lite.py`)
3. **Analyze**: ネットリストを入力し「RUN ANALYSIS」をクリックすると解析結果が表示されます。
   (Enter your netlist and click "RUN ANALYSIS" to view the results)

## ライセンス / License
MIT License
Copyright (c) 2026 zyutama
