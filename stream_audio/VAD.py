import webrtcvad
import numpy as np
import wave
from scipy.io import wavfile


input_file = "/home/kienpt/Documents/Beam/stream_audio/wav/noisy_1588045334.wav"
fs, audio = wavfile.read(input_file)
num_samples, num_channels = audio.shape

win_dur = 0.05
hop_dur = 0.025
vad_win_dur = 0.01

win_size = int(fs*win_dur)
hop_size = int(fs*hop_dur)
vad_win_size = int(fs*vad_win_dur)

vad = webrtcvad.Vad()
vad.set_mode(0)

num_frames = int((num_samples - win_size) / hop_size)

result = np.zeros((num_frames, num_channels), dtype='int16')

for channel in range(num_channels):
    for i in range(num_frames):
        frame = audio[i*hop_size: i*hop_size + win_size, channel]
        res = 0
        for f in range(5):
            begin = f * vad_win_size
            end = (f+1)*vad_win_size
            vad_frame = frame[begin:end].tobytes()
            res += int(vad.is_speech(vad_frame, fs))
        result[i, channel] = round(res/5)

    for i in range(num_frames):
        result[i, :] = round(1.0*sum(result[i, :])/num_channels)

print(result)
np.savetxt('VAD_0.txt',  np.transpose(result), fmt='%d')
