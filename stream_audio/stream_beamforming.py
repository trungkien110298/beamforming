import pyaudio
import webrtcvad
import numpy as np
import wave
import sys
from scipy.io import wavfile
import subprocess


class Speech:
    def __init__(self):
        # Config
        self.chunk = 512
        self.format = pyaudio.paInt16
        self.num_channels = 8
        self.rate = 16000
        self.record_seconds = 3

    def record(self):
        # Start stream
        p = pyaudio.PyAudio()
        stream = p.open(format=self.format,
                        channels=self.num_channels,
                        rate=self.rate,
                        input=True,
                        input_device_index=6,
                        frames_per_buffer=self.chunk)

        print("* recording")
        self.frames = []
        for i in range(0, int(self.rate / self.chunk * self.record_seconds)):
            data = stream.read(self.chunk, exception_on_overflow=False)
            self.frames.append(data)
        print("* done recording")

        framesAll = b''.join(self.frames)
        self.audio = np.fromstring(framesAll, dtype=np.int16)
        chunk_length = int(len(self.audio) / self.num_channels)
        self.audio = np.reshape(self.audio, (chunk_length, self.num_channels))
        self.num_samples = chunk_length

        # Stop stream
        stream.stop_stream()
        stream.close()
        p.terminate()

    # Write audio to file
    def write_to_file(self, file_name, type='raw'):

        # wf = wave.open(WAVE_OUTPUT_FILENAME, 'wb')
        # wf.setnchannels(CHANNELS)
        # wf.setsampwidth(p.get_sample_size(FORMAT))
        # wf.setframeself.rate(self.rate)
        # wf.writeframes(b''.join(frames))
        # wf.close()

        # Write to file
        if type == 'raw':
            wavfile.write(file_name, self.rate, self.audio)
        elif type == 'enhanced':
            wavfile.write(file_name, self.rate, self.audio_enhanced)

    # Calculate VAD
    def cal_VAD(self):
        win_dur = 0.05
        hop_dur = 0.025
        vad_win_dur = 0.01

        win_size = int(self.rate*win_dur)
        hop_size = int(self.rate*hop_dur)
        vad_win_size = int(self.rate*vad_win_dur)

        vad = webrtcvad.Vad()
        vad.set_mode(1)

        sample_per_frame = int(self.rate*10/1000)
        self.num_frames = int((self.num_samples - win_size) / hop_size)

        self.vad = np.zeros(
            (self.num_frames, self.num_channels), dtype='int16')
        for channel in range(self.num_channels):
            for i in range(self.num_frames):
                frame = self.audio[i*hop_size: i*hop_size + win_size, channel]
                res = 0
                for f in range(5):
                    begin = f * vad_win_size
                    end = (f+1)*vad_win_size
                    vad_frame = frame[begin:end].tobytes()
                    res += int(vad.is_speech(vad_frame, self.rate))

                self.vad[i, channel] = round(res/5)

        for i in range(self.num_frames):
            self.vad[i, :] = round(1.0*sum(self.vad[i, :])/self.num_channels)

    # Call C++ module to enhance speech by beamforming
    def beamforming(self):
        # np.set_printoptions(threshold=sys.maxsize)
        audio_str = np.array2string(self.audio, threshold=sys.maxsize)
        audio_str = audio_str.replace('[', '').replace(']', '')

        vad_str = np.array2string(self.vad)
        vad_str = vad_str.replace('[', '').replace(']', '')

        process = subprocess.Popen(
            ["/home/kienpt/Documents/Beam/beamforming_cpp/beamforming"], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        #process.stdin.write(audio_str + "\n"+ vad_str)
        input_str = '\n'.join([str(self.num_samples), str(self.num_channels),
                              str(self.rate), str(self.num_frames), audio_str, vad_str])
        # with open('input.txt', 'w') as f:
        #     f.write('{}'.format(input_str))
        input_bytes = bytes(input_str, 'ascii')
        output_bytes = process.communicate(input=input_bytes)[0]
        output = output_bytes.decode("utf-8")
        output_arr = output.split()

        self.audio_enhanced = np.array(output_arr, dtype=float)
        #self.audio_enhanced = self.audio_enhanced.reshape((self.num_samples, self.num_channels))
        process.stdin.close()


def main():

    s = Speech()
    s.record()
    s.cal_VAD()
    s.beamforming()
    s.write_to_file('python_out.wav', type='enhanced')


if __name__ == "__main__":
    main()
