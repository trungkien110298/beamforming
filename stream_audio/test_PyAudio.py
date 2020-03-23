import pyaudio
import webrtcvad
import numpy as np
import wave
from scipy.io import wavfile


# Config 
CHUNK = 512
FORMAT = pyaudio.paInt16
CHANNELS = 8
RATE = 16000
RECORD_SECONDS = 3
WAVE_OUTPUT_FILENAME = "output.wav"


# Start stream
p = pyaudio.PyAudio()
stream = p.open(format=FORMAT,
                channels=CHANNELS,
                rate=RATE,
                input=True,
                input_device_index= 6,
                frames_per_buffer=CHUNK)

print("* recording")

frames = []
for i in range(0, int(RATE / CHUNK * RECORD_SECONDS)):
    data = stream.read(CHUNK, exception_on_overflow=False)
    frames.append(data)

print(len(frames[0]), len(frames))
print("* done recording")


# Stop stream
stream.stop_stream()
stream.close()
p.terminate()




wf = wave.open(WAVE_OUTPUT_FILENAME, 'wb')
wf.setnchannels(CHANNELS)
wf.setsampwidth(p.get_sample_size(FORMAT))
wf.setframerate(RATE)
wf.writeframes(b''.join(frames))
wf.close()

#Write to file
framesAll = b''.join(frames)
result = np.fromstring(framesAll, dtype=np.int16)
chunk_length = len(result) / CHANNELS
result = np.reshape(result, (chunk_length, CHANNELS))
wavfile.write("wf" + WAVE_OUTPUT_FILENAME,RATE,result)

# Test VAD 

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
