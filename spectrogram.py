# -*- coding: utf-8 -*-
"""
Created on Sat May 21 14:18:37 2016

@author: zhitao wang
@email: zt.wang@hotmail.com
"""

'''
   Module for feature tracking in the dynamic spectrum'''
   
__version__ = '0'

import numpy as np
import matplotlib.pylab as plt
import time
from datetime import datetime
#from scipy import ndimage

class Spectrogram:
    def __init__(self, date = "2012-03-03", data_file = "t18_52_00.00~18_53_00.01_ds_LL.txt", \
            time_file = "t18_52_00.00~18_53_00.01_time_LL_axis.txt", \
            freq_file = "t18_52_00.00~18_53_00.01_freq_LL_axis.txt",\
            file_path = "C:/Users/zhitao/Desktop/pywork/", timerange = None, chanrange = "0~1031", bchan = 0):
        # raw data file and file path
        self.data_file = data_file
        self.time_file = time_file
        self.freq_file = freq_file
        self.file_path = file_path
        # user specified timerange   
        self.timerange = timerange
        self.chanrange = chanrange
        self.date = date
        self.spectrogram = {}
    
    def read_text(self, infile, file_path, dtype = None, delimiter = None):
        ''' Usage: Read data from infile.
            Return: numpy array
        '''
        data = np.loadtxt(file_path+infile, dtype = dtype, delimiter = delimiter)
        return data
    
    def isTimeFormat(self, input):
        ''' Usage: Time format checker
            input: string
            Return: True if input is HHMMSS format, otherwise return False
        '''
        try:
            time.strptime(input, '%H:%M:%S')
            return True
        except ValueError:
            try:
                time.strptime(input, '%H:%M:%S.%f')
                return True
            except ValueError:
                return False    
     
    def read_spectrogram(self, timerange = None, chanrange = None):
        ''' Usage: convert input text to dictionary format.
            timerange: time index or HH:MM:SS, input form: "btime~etime"
            chanrange: absolute channel ID, input form: "0~1031"
            Return: None
        '''
        import re
        
        data = self.read_text(infile = self.data_file, dtype = 'float', file_path = self.file_path)
        time_axis = self.read_text(infile = self.time_file, dtype = 'str', \
            file_path = self.file_path, delimiter = "\n")
        freq_axis = self.read_text(infile = self.freq_file, dtype = 'float', file_path = self.file_path)
        # convert time to seconds
        timesec = np.zeros(len(time_axis), dtype = float)
        # remove "b" character if running in Python 3.x
        for i in range(len(time_axis)):
            time_axis[i] = re.sub("[b']", '', time_axis[i])
            strtime = (self.date + ' ' + time_axis[i])
            dt = datetime.strptime(strtime, '%Y-%m-%d %H:%M:%S.%f')
            timesec[i] = time.mktime(dt.timetuple())
        # unit conversion from Hz to MHz
        freq_axis /= 1e6
        # create CASA format dictionary
        self.spectrogram = {"data": data, "axis_info": \
            {"time_info": {"time_axis": time_axis, "timesec": timesec, "Date": self.date, "unit": "UT (second)"},\
             "freq_info":{ "freq_axis": freq_axis, "unit": "MHz"}}}
        
        self.spectrogram = self.split(timerange = timerange, chanrange = chanrange)
        self.chanrange = chanrange
        # set timerange and chanrange
        self.timerange = self.spectrogram["axis_info"]["time_info"]["time_axis"][0] + \
            '~' + self.spectrogram["axis_info"]["time_info"]["time_axis"][-1]
        
    def plot_spectrogram(self, fig = None, timerange = None, chanrange = None, clim = [1, 100],\
            figsize = None, cmap = 'Greys_r', subplot = 111, aspect = 1):
        ''' Usage: plot object's spectrogram
            fig: plot the current window if fig is specified, else create a new window
            timerange: time index or HH:MM:SS, input form: "btime~etime"
            chanrange: absolute channel ID, input form: "0~1031"
            clim: minimum and maximum value of plot, input form: [min, max]
            figsize: figure size, input form: (xsize, ysize)
            cmap: color table, by default in grey scale
            Return: figure and axis
        '''
        spectrogram = self.split(timerange = timerange, chanrange = chanrange)
        data = spectrogram["data"]
        data = np.flipud(data)
        time_axis = spectrogram["axis_info"]["time_info"]["time_axis"]
        freq_axis = spectrogram["axis_info"]["freq_info"]["freq_axis"]
        freq_unit = spectrogram["axis_info"]["freq_info"]["unit"]
        timerange = [0, len(time_axis)-1]
        freqrange = [freq_axis[0], freq_axis[-1]]
        if fig == None:      
            fig = plt.figure(figsize=figsize)
        else:
            fig = fig
        ax = fig.add_subplot(subplot)
        ax.set_title('Start Time '+time_axis[0])
        plt.xlabel('UT (50 ms/pixel)',{'fontsize':15})
        plt.ylabel('Frequency ('+freq_unit+')',{'fontsize':15})
        plt.imshow(data, extent = timerange + freqrange, cmap = cmap, aspect = aspect, \
            clim = (clim[0], clim[1]))
        return fig, ax
        
    def split(self, timerange = None, chanrange = None):
        ''' Usage: return the spectrogram
            timerange: string format either in HH:MM:SS or pixel time index, joint by "~"
            chanrange: channel range joint by "~"
            Return: splitted data in dictionary format
        '''
        data = self.spectrogram["data"]
        timesec = self.spectrogram["axis_info"]["time_info"]["timesec"]
        time_axis = self.spectrogram["axis_info"]["time_info"]["time_axis"]
        time_unit = self.spectrogram["axis_info"]["time_info"]["unit"]
        freq_axis = self.spectrogram["axis_info"]["freq_info"]["freq_axis"]
        freq_unit = self.spectrogram["axis_info"]["freq_info"]["unit"]
        ny, nx = data.shape
        try:
            timerange = timerange.split('~')
            # convert to pixel time index if it is formated in HH:MM:SS or HH:MM:SS.SS
            for i in range(len(timerange)):
                if self.isTimeFormat(timerange[i]):
                    strtime = (self.date + ' ' + timerange[i])
                    try:
                        dt = datetime.strptime(strtime, '%Y-%m-%d %H:%M:%S')
                    except ValueError:
                        try:                        
                            dt = datetime.strptime(strtime, '%Y-%m-%d %H:%M:%S.%f')
                        except ValueError:
                            raise NameError('Input time format not correct!')      
                    timerange[i] = time.mktime(dt.timetuple())
                else:
                    timerange[i] = timesec[int(timerange[i])]
            btidx, etidx = (np.abs(timesec-timerange[0])).argmin(), (np.abs(timesec-timerange[1])).argmin()
        except (AttributeError, TypeError):
            btidx, etidx = 0, nx
        try:
            chanrange = chanrange.split('~')
            bchan = int(self.chanrange.split('~')[0])
            bchidx, echidx = max(int(chanrange[0])-bchan, 0), min(int(chanrange[1])-bchan, ny-1)
        except (AttributeError, TypeError, ValueError):
            bchidx, echidx = 0, ny
        data = data[bchidx:echidx, btidx:etidx]
        timesec = timesec[btidx:etidx]
        time_axis = time_axis[btidx:etidx]
        freq_axis = freq_axis[bchidx:echidx]
        spectrogram = {"data": data, "axis_info": \
            {"time_info": {"time_axis": time_axis, "timesec": timesec, "Date": self.date, "unit": time_unit},\
                           "freq_info":{ "freq_axis": freq_axis, "unit": freq_unit}}}
        return spectrogram
        
        
    def segmentation(self, threshold):
        img = self.spectrogram["data"]
        mask = (img > threshold).astype(np.float)
        hist, bin_edges = np.histogram(img, bins=60)
        bin_centers = 0.5*(bin_edges[:-1] + bin_edges[1:])
        binary_img = mask > 0.5
        plt.figure(figsize=(11,8))
        plt.subplot(131)
        plt.imshow(img)
        plt.axis('off')
        plt.subplot(132)
        plt.plot(bin_centers, hist, lw=2)
        print(threshold)
        plt.axvline(threshold, color='r', ls='--', lw=2)
        plt.text(0.57, 0.8, 'histogram', fontsize=20, transform = plt.gca().transAxes)
        plt.text(0.45, 0.75, 'threshold = '+ str(threshold)[0:5], fontsize=15, transform = plt.gca().transAxes)
        plt.yticks([])
        plt.subplot(133)     
        plt.imshow(binary_img)
        plt.axis('off')
        plt.subplots_adjust(wspace=0.02, hspace=0.3, top=1, bottom=0.1, left=0, right=1)
        plt.show()
        print(img.max())
        print(binary_img.max())
        
        return mask
        
    def new_spectrum_obj(self, date = "2012-03-03", data_file = "t18_52_00.00~18_53_00.01_ds_LL.txt", \
            time_file = "t18_52_00.00~18_53_00.01_time_LL_axis.txt", \
            freq_file = "t18_52_00.00~18_53_00.01_freq_LL_axis.txt",\
            file_path = "C:/Users/zhitao/Desktop/pywork/", timerange = None, chanrange = None):
        ''' Usage: create a new spectrogram object from the old one.
            timerange: time index or HH:MM:SS, input form: "btime~etime"
            chanrange: absolute channel ID, input form: "0~1031"
            Return: a new spectrogram object with splitted data
        '''
        new_spectrum = Spectrogram()
        new_spectrum.data_file = data_file
        new_spectrum.time_file = time_file
        new_spectrum.freq_file = freq_file
        new_spectrum.date = date
        new_spectrum.spectrogram = self.spectrogram
        new_spectrum.spectrogram = new_spectrum.split(timerange = timerange, chanrange = chanrange)
        new_spectrum.timerange = new_spectrum.spectrogram["axis_info"]["time_info"]["time_axis"][0] + \
            '~' + self.spectrogram["axis_info"]["time_info"]["time_axis"][-1]
        new_spectrum.chanrange = str(self.chanrange.split('~')[0]) + '~' + \
            str(int(self.chanrange.split('~')[0])+ \
            len(new_spectrum.spectrogram["axis_info"]["freq_info"]["freq_axis"]))
        return new_spectrum
        

##################################### code execution/syntax examples below ###################   
######################## initialize an empty Spectrogram object ##############################

#sample = Spectrogram()
#sample = Spectrogram(data = '2012-03-03', timerange = '18:52:00~18:53:00')

######################## read an empty Spectrogram object ####################################

#sample.read_spectrogram(chanrange = '0~1031')
#sample.read_spectrogram(timerange = '18:52:10.02~18:52:19.97', chanrange = '200~800')

########## store the data within specified timerange and chanrange into the sample object ######

#sample.read_spectrogram(chanrange = '200~800')
#threshold = sample.spectrogram["data"].mean()
#mask = sample.segmentation(threshold = threshold)
#img = mask*sample.spectrogram["data"]
#plt.figure(figsize=(11,8))
#plt.imshow(img)
#plt.imshow(mask,cmap=plt.cm.gray)

#sample.plot_spectrogram(value_range = [None, None])


########################### use plot method to view spectrogram ################################

#sample.plot_spectrogram(value_range = [0.3, 1.5], timerange = '18:52:10.02~18:52:19.97', chanrange = '200~800')
#sample.plot_spectrogram(value_range = [0.3, 1.5], timerange = '18:52:10.02~18:52:19.97', chanrange = '200~500')
#sample.plot_spectrogram(value_range = [0.3, 1.5], timerange = '18:52:10.02~18:52:19.97', chanrange = '500~800')

################### you can also create another object from the original object ################

#sample2 = sample.new_spectrum_obj()
#f, ax = sample2.plot_spectrogram(timerange = '18:52:10.02~18:52:19.97', chanrange = '200~800')
    
    