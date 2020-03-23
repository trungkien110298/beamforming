import webrtcvad
import numpy as np
import wave
from scipy.io import wavfile



input_file = "../data/speaker0250-0020.wav"
fs, audio = wavfile.read(input_file)
num_samples, num_channels = audio.shape

win_size = 0.01
samples_per_frame = int(fs*win_size)

# init VAD
vad = webrtcvad.Vad()
vad.set_mode(1)


# calculate VAD
result = np.zeros((num_samples,num_channels), dtype='int16')
for channel in range(num_channels):
    begin = 0
    print("Channel", channel)
    while begin < num_samples:
        end = min(begin + samples_per_frame, num_samples)
        frame = audio[begin:end, channel]
        if frame.shape[0] < samples_per_frame:
            remain_len = samples_per_frame-frame.shape[0]
            frame = np.concatenate((frame, np.zeros(remain_len, dtype='int16')), axis=0)
        
        frame = frame.tobytes()
        t = vad.is_speech(frame, fs)
       
        result[begin:end, channel].fill(int(t)) 
        begin = begin + samples_per_frame

np.savetxt('VAD.txt',  np.transpose(result), fmt= '%d')
        

