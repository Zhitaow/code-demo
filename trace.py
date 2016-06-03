# -*- coding: utf-8 -*-
"""
Created on Wed May 25 15:42:11 2016

@author: zhitao wang
@email: zt.wang@hotmail.com
"""

'''
   Description: modules for purpose of analyzing fine structures in the dynamic spectrum.
   function: 1. automatic loop tracing on intermediate drift bursts
             2. statistical analysis on traced bursts
   '''
   
__version__ = '0'

import numpy as np
import scipy.ndimage as ndimage
import matplotlib.pylab as plt
from spectrogram import Spectrogram

class Trace:
    
    def __init__(self, para = [2, 800, 100, 1500, 0., 1000], outfile = 'trace.txt'):
        self.window_len = int(para[0])  # smoothing window length (pixel unit)
        self.rmin = int(para[1])        # minimum curvature radius
        self.lmin = int(para[2])        # minimum of traced loop length
        self.nstruc = int(para[3])      # maximal number of structures for analysis
        self.qmed = para[4]             # threshold ratio
        self.nloopmax = para[5]         # maximal number of loop
        self.outfile = outfile          # None: no outputfile is generated
        self.trace_info = {}            # the coordinates of traced structures 
        
    def nanargmax(self,a):
        ''' Usage: find the maxiumum of array, exclude nan value
             Return: the maxium value and the indice as tuple
        '''
        idx = np.argmax(a, axis=None)
        multi_idx = np.unravel_index(idx, a.shape)
        if np.isnan(a[multi_idx]):
            print("Ignore NANs in maxtrice. ")
            nan_count = np.sum(np.isnan(a))
            # In numpy < 1.8 use idx = np.argsort(a, axis=None)[-nan_count-1]
            idx = np.argpartition(a, -nan_count-1, axis=None)[-nan_count-1]
            multi_idx = np.unravel_index(idx, a.shape)
        max_value = a[multi_idx]
        return max_value, multi_idx
        
    def trace(self, spectrogram, mask = None):
        ''' Usage: tracing algorithm for stripe pattern in the spectrogram
            Return: the coordinates of traced structures 
            
        '''
        window_len = self.window_len
        rmin = self.rmin
        lmin = self.lmin
        nstruc = self.nstruc
        qmed = self.qmed
        nloopmax = self.nloopmax
        outfile = self.outfile
        data = spectrogram["data"]
        dim = data.shape
        nx = dim[0]
        ny = dim[1]
        timesec = spectrogram["axis_info"]["time_info"]["timesec"]
        time_axis = spectrogram["axis_info"]["time_info"]["time_axis"]
        freq_axis = spectrogram["axis_info"]["freq_info"]["freq_axis"]
        ########################  apply mask to data  #########################
        try:
            mx, my = mask.shape
            if mx != nx or my != ny:
                raise ValueError("Mask's dimension not matched with spectrogram!")
        except (AttributeError, TypeError):
            # if mask not specified, do not apply the mask
            mask = np.ones(data.shape)
        data *= mask
        ############################ define some constants ####################
        reso = 1
        step = 1
        ngap = 3
        npmax = 1000
        window_len_ = window_len + 2
        nlen = rmin
        na = 180
        alpha = np.pi*np.arange(0,na)/float(na)
        ndir = 2 
        nb = 30
        s_loop = step*np.arange(nlen)
        s0_loop = step*(np.arange(nlen)-nlen/2)
        wid	= max(window_len_/2-1, 1)
        # image base level
        fluxmin = data.min()
        fluxmax = data.max()
        fluxmed = np.median(data)
        base = np.ones(dim)+fluxmed*qmed
        print("Min/Max flux in spectrogram = " + str(fluxmin) + ',' + str(fluxmax))
        data1 = np.maximum(data, base)
        
        # Highpass filter
        if window_len <= 2:
            data2 = data1 - ndimage.gaussian_filter(data1, sigma=[window_len_]*2)
        else:
            data2 = ndimage.gaussian_filter(data1, sigma=[window_len]*2) \
                - ndimage.gaussian_filter(data1, sigma=[window_len_]*2)
        # remove boundary zone due to smoothing effect
        data2[0:window_len_,0:-1] = 0
        data2[nx-window_len_:nx,0:-1] = 0
        data2[0:-1,0:window_len_] = 0
        data2[0:-1,ny-window_len_:ny] = 0
        
        dim = data2.shape
        nx = dim[0]
        ny = dim[1]
        residue = np.maximum(data2, np.zeros(data2.shape))
        threshold = np.median(residue)
        data0 = np.copy(residue)
        print("The threshold is: "+str(threshold))

        iloop = 0 # tracing loop ID
        ilast = -1 
        iloop_nstruc = np.zeros(nstruc)

        # trace starting point at the maximum-value pixel of the image
        for istruc in range(nstruc):
            skip = False
            zstart, (istart, jstart) = self.nanargmax(residue)
            if zstart <= 0:
                break
            if istruc % 100 == 0:
                # estimate residue area in the map
                npix = np.count_nonzero(residue)
                qpix = float(npix)/(nx*ny)*100
                print("Struc #" + str(istruc) + " ID #" + str(iloop) + " Flux = " \
                    + str(zstart)[0:5] + " Residue area = " + str(qpix)[0:5] +" %")
                if iloop == ilast:
                    break
                ilast = iloop
            # Tracing stepwise
            ip = 0
            # loop body for bi-direction
            for idir in range(ndir):
                [xl, yl, zl, al, ir] = [np.zeros(npmax + 1) for array in range(5)]
                if idir == 0:
                    dir_sign = 1
                elif idir == 1:
                    dir_sign = -1
                # find initial direction
                xl[0], yl[0], zl[0] = istart, jstart, zstart
                flux_max = 0.
                # loop body for initial angular direction
                for ia in range(int(na)):
                    x_ = xl[0] + s0_loop*np.cos(alpha[ia])
                    y_ = yl[0] + s0_loop*np.sin(alpha[ia])
                    ix = (x_ + 0.5).astype(int)
                    iy = (y_ + 0.5).astype(int)
                    # confine tracing coordinates within the data
                    idx = np.logical_and(ix >= 0, ix < nx)
                    ix, iy = ix[idx], iy[idx]
                    idy = np.logical_and(iy >= 0, iy < ny)
                    ix, iy = ix[idy], iy[idy]
                    flux_ = np.maximum(residue[ix,iy],0)
                    flux = np.sum(flux_)/float(nlen)
                    if flux > flux_max:
                        flux_max = flux
                        al[0] = alpha[ia]
                
                # loop body for maximum loop length
                for ip in np.arange(0,npmax):
                    if ip == 0:
                        # large range of curv radius
                        ib1, ib2 = 0, nb-1
                    else:
                        # small range of curv radius
                        ib1, ib2 = int(np.maximum(ir[ip]-1, 0)), int(np.minimum(ir[ip]+1, nb-1))
                    # angle at curvature center
                    beta0 = al[ip] + np.pi/2.
                    xcen = xl[ip] + rmin*np.cos(beta0)
                    ycen = yl[ip] + rmin*np.sin(beta0)
                    flux_max = 0.
                    # loop for trying different curvature
                    for ib in np.arange(ib1, ib2+1):
                        rad_i = rmin/(-1. + 2.*float(ib)/float(nb-1))
                        xcen_i = xl[ip] + (xcen - xl[ip])*(rad_i/rmin)          # position of curvature center
                        ycen_i = yl[ip] + (ycen - yl[ip])*(rad_i/rmin)
                        beta_i = beta0 + dir_sign*s_loop/rad_i                  # angle at curvature center
                        x_ = xcen_i - rad_i*np.cos(beta_i)                      # x-coordinate of curved segment
                        y_ = ycen_i - rad_i*np.sin(beta_i)                      # 
                        ix = (x_ + 0.5).astype(int)                             # x index of curved segment
                        iy = (y_ + 0.5).astype(int)
                        # confine tracing coordinates within the data
                        idx = np.logical_and(ix >= 0, ix < nx)
                        ix, iy = ix[idx], iy[idx]
                        idy = np.logical_and(iy >= 0, iy < ny)
                        ix, iy = ix[idy], iy[idy]
                        flux_ = np.maximum(residue[ix,iy],0)
                        flux = np.sum(flux_)/float(nlen)
                                            
                        if flux > flux_max:
                            flux_max = flux
                            al[ip+1] = al[ip] + dir_sign*(step/rad_i)
                            ir[ip+1] = ib                                       # curvature radius index
                            al_mid = (al[ip] + al[ip+1])/2.
                            xl[ip+1], yl[ip+1] = xl[ip] + step*np.cos(al_mid + np.pi*idir), yl[ip] + step*np.sin(al_mid + np.pi*idir)                                                   
                            ix_ip, iy_ip = np.clip(int(xl[ip+1] + 0.5), 0, nx-1), np.clip(int(yl[ip+1] + 0.5), 0, ny-1)
                            zl[ip+1] = residue[ix_ip, iy_ip]

                    iz1 = np.maximum(ip+1-ngap, 0)
                    if zl[iz1:ip+2].max() <= 0:
                        ip = np.maximum(iz1-1, 0)
                        break # end segm
                # re-order the trace x-y coordinates
                if idir == 0:
                    # reverse the order in one tracing direction
                    xloop, yloop, zloop = xl[0:ip+1][::-1], yl[0:ip+1][::-1], zl[0:ip+1][::-1]
                elif idir ==1 and ip >=1:
                    # adjoin the other half with opposive direction
                    xloop, yloop, zloop = np.concatenate((xloop, xl[1:ip+1])), np.concatenate((yloop, yl[1:ip+1])), np.concatenate((zloop, zl[1:ip+1]))
            # end for idir =0, 1
                    
            boolarr = np.logical_and(xloop != 0, yloop != 0)
            ind = np.where(boolarr == True)
            nind = len(ind[0])
            looplen = 0
            
            # skip tracking if it's below threshold or exceed the maxmium number of loops
            if nind <= 1:
                skip = True
            else:
                xloop, yloop, zloop = xloop[ind], yloop[ind], zloop[ind]
                if np.median(zloop) <= threshold:
                    skip = True
                else:
                    skip = False
                if iloop >= nloopmax:
                    break                                                       # end structure tracing
   
            if not skip:
                NP = len(xloop)
                s = np.zeros(NP)
                looplen = 0
                if NP >= 2:
                    for ip in np.arange(1, NP):
                        s[ip] = s[ip-1] + np.sqrt((xloop[ip] - xloop[ip-1])**2+(yloop[ip] - yloop[ip-1])**2)
                looplen = s[NP-1]                                               # number of pixels for full loop length
                ns = np.maximum(int(looplen),3)

                # store loop coordinates
                if looplen >= lmin:
                    nn = int(ns/reso +0.5)
                    # save tracing coordinates
                    if iloop == 0:
                        trace_id, trace_x, trace_y, trace_flux = np.tile(iloop, len(xloop)), xloop, yloop, zloop
                    else:
                        trace_id, trace_x, trace_y, trace_flux = np.concatenate((trace_id, np.tile(iloop, len(xloop)))), np.concatenate((trace_x, xloop)),\
                            np.concatenate((trace_y, yloop)), np.concatenate((trace_flux, zloop))
                    # export to outfile
                    if outfile != None:
                        if iloop == 0:
                            write_method = 'w'
                        else:
                            write_method = 'a'
                        
                        with open(outfile, write_method) as file:
                            for ip in range(nn):
                                line = str(iloop) + ', ' + str(xloop[ip]) + ', ' + str(yloop[ip]) + ', ' + str(zloop[ip])
                                file.write(line+'\n')
                        file.close()
                    
                    iloop_nstruc[istruc] = iloop
                    iloop += 1
                    
            # erase loop in residue spectrogram
            i3 = np.maximum(istart-wid, 0)
            i4 = np.minimum(istart+wid+1, nx)
            j3 = np.maximum(jstart-wid, 0)
            j4 = np.minimum(jstart+wid+1, ny)
            residue[i3:i4,j3:j4] = 0
            nn = len(xloop)
            for i in range(nn):
                i0 = int(np.clip(xloop[i], 0, nx-1))
                i3 = int(np.maximum(i0-wid,0))
                i4 = int(np.minimum(i0+wid+1, nx))

                j0 = int(np.clip(yloop[i], 0, ny-1))
                j3 = int(np.maximum(j0-wid,0))
                j4 = int(np.minimum(j0+wid+1, ny))
                residue[i3:i4, j3:j4] = 0
            # direction
        # structure
        spectrogram = {"data": data0, "axis_info": \
            {"time_info": {"time_axis": time_axis, "timesec": timesec, "unit": "UT (second)"},\
             "freq_info":{ "freq_axis": freq_axis, "unit": "MHz"}}}
            
        trace_info = {"ID":trace_id , "indx": trace_x, "indy": trace_y, "flux": trace_flux, "spectrogram": spectrogram}
        return trace_info
        
        
    def plot_trace(self, spectrogram, trace_info, clim = [0.5, 1.5], color = 'r'):
        ''' plot traced coordinates in the dynamic spectrum
        '''
        
        base0 = Spectrogram()
        base0.spectrogram = spectrogram
        fig0, ax0 = base0.plot_spectrogram(figsize = [12,8], subplot = 211, clim = clim, cmap = 'Greys_r', aspect = 0.2)
        base1 = Spectrogram()
        base1.spectrogram = spectrogram
        fig1, ax1 = base1.plot_spectrogram(fig = fig0, subplot = 212, clim = clim, cmap = 'Greys_r', aspect = 0.2)
        trace_info = self.idx2map(spectrogram, trace_info)
        idplt = trace_info["ID"]
        xplt= trace_info["indy"]
        yplt = trace_info["frequency"]
        # overplot traced drift bursts
        for i in np.arange(idplt.min(), idplt.max()+1):
            idx = np.where(idplt == i)[0]
            if len(idx) > 1:
                plt.plot(xplt[idx],yplt[idx], color = color)
        ylim = ax0.get_ylim()
        ax1.set_ylim(ylim)
        plt.show()
        
        
    def idx2map(self, spectrogram, trace_info):
        ''' map out the time and frequency according to the traced coorindates (indices)
        '''
        timesec = spectrogram["axis_info"]["time_info"]["timesec"]
        time_axis = spectrogram["axis_info"]["time_info"]["time_axis"]
        freq_axis = spectrogram["axis_info"]["freq_info"]["freq_axis"]
        trace_id = trace_info["ID"]
        indx = trace_info["indx"]
        indy = trace_info["indy"]
        flux = trace_info["flux"]
        
        [trace_timesec, trace_freq] = [np.zeros(len(trace_id)) for array in range(2)]
        
        for i in range(len(trace_id)):
            if i == 0:
                trace_time = time_axis[int(indy[i])]
            else:
                trace_time = [trace_time, time_axis[int(indy[i])]]
                
            trace_timesec[i] = timesec[int(indy[i])]
            trace_freq[i] = freq_axis[int(indx[i])]

        trace_info = {'ID':trace_id, "indx": indx, "indy": indy, 'time': trace_time, 'timesec': trace_timesec, 'frequency': trace_freq, 'flux': flux}
        return trace_info
        
###################### code execution/syntax examples below ###################
###########################  read spectrogram #################################
#sample = Spectrogram()
#sample.read_spectrogram(timerange = '18:52:10.02~18:52:19.97', chanrange = '200~800')
#sample.read_spectrogram()
###########################  plot spectrogram #################################
#sample.plot_spectrogram()
#data = sample.spectrogram
###########################  mask spectrogram #################################
#threshold = sample.spectrogram["data"].mean()
#mask = sample.segmentation(threshold = threshold)
##########################  trace spectrogram #################################
#t=Trace()
#trace_info = t.trace(data, mask = None)
#t.plot_trace(data, trace_info)

    