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
1. Download the files above to your directory.
2. Import raw data.

  In[1]: from trace import Trace

  In[2]: from spectrogram import Spectrogram

  In[3]: file_path = "YOUR_PATH"

  In[4]: sample = Spectrogram(file_path = file_path)

  In[5]: sample.read_spectrogram()
  
3. Apply tracing over specific range in time and frequency.
  In[6]: t=Trace()

  In[7]: trace_info = t.trace(data, mask = None)
  
4. Visualize traced structures.
  In[8]: t.plot_trace(data, trace_info)

