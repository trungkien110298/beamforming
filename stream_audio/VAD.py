import webrtcvad
import numpy as np
import wave
from scipy.io import wavfile



WAVE_INPUT_FILENAME = "../data/speaker0250-0020.wav"
fs, audio = wavfile.read(WAVE_INPUT_FILENAME)
print(audio.shape)


vad = webrtcvad.Vad()
vad.set_mode(1)

begin = 0
sample_per_frame = int(RATE*10/1000)

while True:
    end = begin + sample_per_frame
    if end > result.shape[0]:
        break
    test_frame = (result[begin:end,0]).tobytes()
    print(vad.is_speech(test_frame, RATE))
    begin = begin + sample_per_frame

