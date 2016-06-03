# -*- coding: utf-8 -*-
"""
Created on Thu Jun 02 23:35:05 2016

@author: Zhitao
"""

from trace import Trace
from spectrogram import Spectrogram
########################## Import raw data. ###################################
file_path = "C:/Users/Zhitao/Desktop/pywork/"           # modify YOUR_PATH here
sample = Spectrogram(file_path = file_path)
sample.read_spectrogram()
data = sample.spectrogram
########### Apply tracing over specific range in time and frequency. ##########
t = Trace()
trace_info = t.trace(data, mask = None)
# Visualize traced structures.
t.plot_trace(data, trace_info)
