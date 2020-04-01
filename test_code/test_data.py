from scipy.io import wavfile

__, samples = wavfile.read("/home/kienpt/Documents/Beam/data/speaker0250-0020.wav")
print(samples)
