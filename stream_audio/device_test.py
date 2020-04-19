import pyaudio


pa = pyaudio.PyAudio()
print(pa.get_default_input_device_info())
print(pa.get_device_count())
for i in range (pa.get_device_count()):
    print(pa.get_device_info_by_index(i))
    print("\n")
