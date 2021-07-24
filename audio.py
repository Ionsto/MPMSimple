import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as signal
from scipy.io.wavfile import write

def write_file(data_raw,Fs=44100):
    step_original = 1/Fs
    step_new = 1/44100.0
    t_end = len(data_raw) / Fs

    t_original = np.arange(0,t_end,step_original)
    t_new = np.arange(0,t_end,step_new)
    data = np.interp(t_new,t_original,data_raw)
    plt.figure()
    plt.plot(t_original,data_raw)
    plt.plot(t_new,data)
    scaled = np.int16(data/np.max(np.abs(data)) * 32767)
    write('test.wav', 44100, scaled)

data = np.loadtxt("audio_samples.txt")



dt = data[0]
data = data[1:]
Fs = 1/dt

sos = signal.butter(10, 50, 'hp', fs=Fs, output='sos')

filtered = signal.sosfilt(sos,data)
write_file(filtered,Fs)
