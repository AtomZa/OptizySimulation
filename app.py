import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.signal import find_peaks

import glob
import  os

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

# Setting
st.set_option('deprecation.showPyplotGlobalUse', False)

# Navigation
rad = st.sidebar.radio("Menu",["Preprocessing","Analysis","Result"])

# Upload Files
uploaded_files = st.sidebar.file_uploader("Choose a CSV file", accept_multiple_files=True)

# Parameter 1
st.sidebar.subheader("Smooth & Interpolate")
ns = st.sidebar.number_input("Smooth Factor", min_value = 0, value = 100, step = 10)
ir = st.sidebar.number_input("Interpolate Interval", min_value = 0.000, value = 0.001, step = 0.001, format="%.3f")

# Parameter 2
st.sidebar.subheader("Time Interval")
c, d, e = st.columns(3)
starttime = st.sidebar.number_input("Start Time", min_value = 0, value = 1, step = 1)
timeinterval = st.sidebar.number_input("Time Interval", min_value = 0, value = 5, step = 1)
endtime = st.sidebar.number_input("End Time", min_value = 0, value = 26, step = 1)
timeaxis = np.arange(starttime, endtime+timeinterval, timeinterval)

# Parameter 3
st.sidebar.subheader("Wavelength")
minwavelength = st.sidebar.number_input("Min Wavelength", min_value = 0, value = 200, step = 10)
maxwavelength = st.sidebar.number_input("Max Wavelength", min_value = 0, value = 1000, step = 10)

if rad == "Preprocessing":
    st.title("Preprocessing")
    # Raw Data
    path = "testdata"
    filelist = glob.glob(os.path.join(path, '*.csv'))
    for idx, uploaded_file in enumerate(uploaded_files):
        if idx < starttime: pass
        elif idx > endtime: break
        elif (idx-1) % timeinterval != 0: pass
        else:
            df = pd.read_csv(uploaded_file, header=None)
            df.columns = ["wavelength", "intensity"]
            df["intensity"].plot()
    st.pyplot()
        
if rad == "Analysis": 
    st.title("Analysis")
    # Spectrum 
    smootht1y_max = np.zeros(len(timeaxis))
    index = np.zeros(len(timeaxis))
    b = 0

    for idx, uploaded_file in enumerate(uploaded_files):
        if idx < starttime: pass
        elif idx > endtime: break
        elif (idx-1) % timeinterval != 0: pass
        else:
            df = pd.read_csv(uploaded_file, header=None)
            df.columns = ["wavelength", "intensity"]
            t1x = df["wavelength"]
            t1y = df["intensity"]

            filtert1y = df.loc[(df['wavelength'] < maxwavelength) & (df['wavelength'] > minwavelength)]
            filtert1y = filtert1y['intensity'].to_numpy().flatten()
            filtert1x = np.linspace(minwavelength, maxwavelength, max(filtert1y.shape))

            smootht1y = smooth(filtert1y,ns)

            interpt1y = interp1d(filtert1x,smootht1y)
            interpt1x = np.linspace(minwavelength, maxwavelength, round(max(filtert1y.shape)*ir*100))

            smootht1y_max[b] = max(interpt1y(interpt1x))
            index[b] = np.argmax(interpt1y(interpt1x))
            b = b+1
            plt.plot(interpt1y(interpt1x))

            peaks, _ = find_peaks(interpt1y(interpt1x), distance=1)
            np.diff(peaks)
            plt.plot(peaks, interpt1y(interpt1x)[peaks], "x", color="red")    
    st.pyplot()

if rad == "Result":
    st.balloons()
    st.title("Result")
    smootht1y_max = np.zeros(len(timeaxis))
    index = np.zeros(len(timeaxis))
    b = 0

    for idx, uploaded_file in enumerate(uploaded_files):
        if idx < starttime: pass
        elif idx > endtime: break
        elif (idx-1) % timeinterval != 0: pass
        else:
            df = pd.read_csv(uploaded_file, header=None)
            df.columns = ["wavelength", "intensity"]
            t1x = df["wavelength"]
            t1y = df["intensity"]

            filtert1y = df.loc[(df['wavelength'] < maxwavelength) & (df['wavelength'] > minwavelength)]
            filtert1y = filtert1y['intensity'].to_numpy().flatten()
            filtert1x = np.linspace(minwavelength, maxwavelength, max(filtert1y.shape))

            smootht1y = smooth(filtert1y,ns)

            interpt1y = interp1d(filtert1x,smootht1y)
            interpt1x = np.linspace(minwavelength, maxwavelength, round(max(filtert1y.shape)*ir*100))

            smootht1y_max[b] = max(interpt1y(interpt1x))
            index[b] = np.argmax(interpt1y(interpt1x))
            b = b+1

            peaks, _ = find_peaks(interpt1y(interpt1x), distance=1)
            np.diff(peaks)

    plt.subplot(1,2,1)
    plt.plot(timeaxis, index, '.')
    plt.subplot(1,2,2)
    plt.plot(timeaxis, smootht1y_max, '.')
    st.pyplot()

