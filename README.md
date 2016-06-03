# code-demo
### This repo contains demo-codes of auto tracing algorithm.

#### raw input data including:
- [x] 1-miniute spectrogram obtained by VLA: 18_52_00.00~18_53_00.01_ds_LL.txt
- [x] x-axis (time) info: t18_52_00.00~18_53_00.01_time_LL_axis.txt
- [x] y-axis (frequency) info: t18_52_00.00~18_53_00.01_freq_LL_axis.txt

#### python program including:
- [x] reading and plotting spectrogram: spectrogram.py
- [x] tracing algorithm: trace.py

### Example: pattern discovery of solar fiber bursts in the radio spectrogram
1. Download the files above into your directory.
2. Import raw data
  file_path = "YOUR_PATH" 
  <li>sample = Spectrogram(file_path = file_path)</li>
  <li>sample.read_spectrogram()</li>
3. Apply tracing over specific range in time and frequency.
  <li>t=Trace()</li>
  <li>trace_info = t.trace(data, mask = None)</li>
4. Visualize traced structures
  <li>t.plot_trace(data, trace_info)</li>
