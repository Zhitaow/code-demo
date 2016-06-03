# code-demo
## This repo contains demo-codes of auto tracing algorithm.

### Description:
Perform graphical texture tracing for use of further data analysis.

### Tracing Result:
![alt tag](https://github.com/Zhitaow/code-demo/blob/master/figure_1.png)

#### raw input data including:
- [x] 1-miniute spectrogram obtained by VLA: "18_52_00.00~18_53_00.01_ds_LL.txt"
- [x] x-axis (time) info: "t18_52_00.00~18_53_00.01_time_LL_axis.txt"
- [x] y-axis (frequency) info: "t18_52_00.00~18_53_00.01_freq_LL_axis.txt"

#### python program including:
- [x] reading and plotting spectrogram: "spectrogram.py"
- [x] tracing algorithm: "trace.py"

### Example: pattern discovery of solar fiber bursts in the radio spectrogram
1. Download the files above to your directory.
2. In "testme.py", specify "file_path" variable corresponding to your directory
3. Run wrap-up script: "testme.py"
