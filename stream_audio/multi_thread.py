import pyaudio
import webrtcvad
import numpy as np
import wave
import sys
from scipy.io import wavfile
import subprocess
import pyloudnorm as pyln
from queue import Queue
from threading import Thread
import time

CHUNK = 512
SAMPLE_RATE = 16000
NUM_CHANNEL_INPUT = 8
NUM_CHANNEL_OUTPUT = 1
BATCH_DURATION = 3
FORMAT = pyaudio.paInt16
INPUT_DEVICE_INDEX = 6


class Record(Thread):
    """This is the worker thread function.
    It record audio form microphone and push to record queue.
    These daemon threads go into an infinite loop,
    and only exit when the main thread ends.
    """

    def __init__(self, format=pyaudio.paInt16, num_channels=8, rate=16000, chunk=512, input_device_index=6, batch_dur=3, record_queue=None, raw_record_queue=None):
        Thread.__init__(self)
        self.format = format
        self.num_channels = num_channels
        self.rate = rate
        self.chunk = chunk
        self.input_device_index = input_device_index
        self.batch_dur = batch_dur
        self.record_queue = record_queue
        self.raw_record_queue = raw_record_queue

    def run(self):
        p = pyaudio.PyAudio()
        stream = p.open(format=self.format,
                        channels=self.num_channels,
                        rate=self.rate,
                        input=True,
                        input_device_index=self.input_device_index,
                        frames_per_buffer=self.chunk)

        sample_per_frame = int(self.batch_dur*self.rate/self.chunk)

        for n in range(3):
            print('record {}'.format(n))
            frames = []
            for i in range(sample_per_frame):
                data = stream.read(self.chunk, exception_on_overflow=False)
                frames.append(data)

            framesAll = b''.join(frames)

            batch_audio = np.fromstring(framesAll, dtype=np.int16)
            chunk_length = int(len(batch_audio) / self.num_channels)
            batch_audio = np.reshape(
                batch_audio, (chunk_length, self.num_channels))
            self.record_queue.put(batch_audio)
            self.raw_record_queue.put(batch_audio)
            # Stop stream

        stream.stop_stream()
        stream.close()
        p.terminate()


class BeamForming(Thread):
    """This is the worker thread function.
    It enhances noisy audio form record queue and push to enhanced queue.
    These daemon threads go into an infinite loop,
    and only exit when the main thread ends.
    """

    def __init__(self, format=pyaudio.paInt16, num_channels=8, rate=16000, chunk=512, input_device_index=6, batch_dur=3, record_queue=None, enhanced_queue=None,  vad_queue=None):
        Thread.__init__(self)
        self.format = format
        self.num_channels = num_channels
        self.rate = rate
        self.chunk = chunk
        self.input_device_index = input_device_index
        self.batch_dur = batch_dur
        self.record_queue = record_queue
        self.enhanced_queue = enhanced_queue
        self. vad_queue = vad_queue

    def run(self):

        # Do not change these configs
        win_dur = 0.05
        hop_dur = 0.025
        vad_win_dur = 0.01

        win_size = int(self.rate*win_dur)
        hop_size = int(self.rate*hop_dur)
        vad_win_size = int(self.rate*vad_win_dur)

        vad = webrtcvad.Vad()
        vad.set_mode(2)

        for n in range(3):
            # Calculate VAD

            # sample_per_frame = int(self.rate*10/1000)

            while True:
                if not self.record_queue.empty():
                    break
            print('beamforming {}'.format(n))
            batch_audio = self.record_queue.get()
            num_samples = batch_audio.shape[0]
            num_frames = int((num_samples - win_size) / hop_size)

            vad_array = np.zeros(
                (num_frames, self.num_channels), dtype='int16')

            for channel in range(self.num_channels):
                for i in range(num_frames):
                    frame = batch_audio[i*hop_size: i *
                                        hop_size + win_size, channel]
                    res = 0
                    for f in range(5):
                        begin = f * vad_win_size
                        end = (f+1)*vad_win_size
                        vad_frame = frame[begin:end].tobytes()
                        res += int(vad.is_speech(vad_frame, self.rate))

                    vad_array[i, channel] = round(res/5)

            for i in range(num_frames):
                vad_array[i, :] = round(
                    1.0*sum(vad_array[i, :])/self.num_channels)
            self.vad_queue.put(vad_array)
            # Beamforming

            # np.set_printoptions(threshold=sys.maxsize)
            batch_audio = 1.0*batch_audio/32767
            audio_str = np.array2string(batch_audio, threshold=sys.maxsize)
            audio_str = audio_str.replace('[', '').replace(']', '')

            vad_str = np.array2string(vad_array)
            vad_str = vad_str.replace('[', '').replace(']', '')

            process = subprocess.Popen(
                ["/home/kienpt/Documents/Beam/beamforming_cpp/beamforming"], stdin=subprocess.PIPE, stdout=subprocess.PIPE)

            input_str = '\n'.join([str(num_samples), str(self.num_channels),
                                   str(self.rate), str(num_frames), str(n+1), audio_str, vad_str])

            input_bytes = bytes(input_str, 'ascii')
            output_bytes = process.communicate(input=input_bytes)[0]
            output = output_bytes.decode("utf-8")
            output_arr = output.split()

            audio_enhanced = np.array(output_arr, dtype=float)
            audio_enhanced = np.array(audio_enhanced*32767, dtype='int16')
            process.stdin.close()
            self.enhanced_queue.put(audio_enhanced)


class play_thread(Thread):
    """This is the worker thread function.
    It gets enhanced audio form enhanced queue and sends to virtual microphone.  
    These daemon threads go into an infinite loop, 
    and only exit when the main thread ends.
    """


def main():

    record_queue = Queue()
    raw_record_queue = Queue()
    enhanced_queue = Queue()
    vad_queue = Queue()

    record = Record(record_queue=record_queue,
                    raw_record_queue=raw_record_queue)

    enhance = BeamForming(record_queue=record_queue,
                          enhanced_queue=enhanced_queue, vad_queue=vad_queue)

    record.start()
    enhance.start()

    record.join()
    enhance.join()

    vad = vad_queue.get()
    while not vad_queue.empty():
        batch_vad = vad_queue.get()
        vad = np.append(vad, batch_vad, axis=0)
    #np.savetxt('VAD.txt',  np.transpose(vad), fmt='%d')

    # Write raw audio
    timestamp = str(int(time.time()))
    if True:
        audio = raw_record_queue.get()
        while not raw_record_queue.empty():
            batch_audio = raw_record_queue.get()
            audio = np.append(audio, batch_audio, axis=0)

        norm_audio = pyln.normalize.peak(audio, -1.0)
        file_name = 'noisy_' + timestamp + '.wav'
        wavfile.write('wav/' + file_name, SAMPLE_RATE, norm_audio)

    if True:
        audio = enhanced_queue.get()
        while not enhanced_queue.empty():
            batch_audio = enhanced_queue.get()
            audio = np.append(audio, batch_audio, axis=0)

        norm_audio = pyln.normalize.peak(audio, -1.0)
        file_name = 'enhanced_' + timestamp + '.wav'
        wavfile.write('wav/' + file_name, SAMPLE_RATE, norm_audio)


if __name__ == "__main__":
    main()
