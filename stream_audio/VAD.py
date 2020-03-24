import webrtcvad
import numpy as np
import wave
from scipy.io import wavfile



input_file = "../data/speaker0250-0020_0.5m.wav"
fs, audio = wavfile.read(input_file)
num_samples, num_channels = audio.shape

win_dur = 0.05
hop_dur = 0.025 
vad_win_dur = 0.01

win_size = int(fs*win_dur)
hop_size = int(fs*hop_dur)
vad_win_size = int(fs*vad_win_dur)
# init VAD
vad = webrtcvad.Vad()
vad.set_mode(1)


# calculate VAD
num_frames = int((num_samples - win_size) / hop_size);
print(num_frames)

result = np.zeros((num_frames, num_channels), dtype='int16')
for channel in range(num_channels):
    begin = 0
    count = 0
    for i in range(num_frames):
        frame = audio[i*hop_size: i*hop_size + win_size, channel]        
        res = 0
        for f in range(5):    
            vad_frame = frame[f*vad_win_size: (f+1)*vad_win_size].tobytes()
            res += int(vad.is_speech(vad_frame, fs))
       
        result[i, channel] = round(res/5) 

for i in range(num_frames):
    result[i, :] = round(1.0*sum(result[i,:])/num_channels)        
print(result.shape)
np.savetxt('VAD.txt',  np.transpose(result), fmt= '%d')
        

