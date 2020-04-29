#!/bin/bash

# This script will create a virtual microphone for PulseAudio to use and set it as the default device.

# Load the "module-pipe-source" module to read audio data from a FIFO special file.
echo "Creating virtual microphone."
pactl load-module module-pipe-source source_name=virtmic file=/home/kienpt/audioFiles/virtmic format=s16le rate=16000 channels=1

# Set the virtmic as the default source device.
echo "Set the virtual microphone as the default device."
pactl set-default-source virtmic

# Create a file that will set the default source device to virtmic for all
#PulseAudio client applications.
echo "default-source = virtmic" > /home/kienpt/.config/pulse/client.conf

# Write the audio file to the named pipe virtmic. This will block until the named pipe is read.
echo "Writing audio file to virtual microphone."

for ((i=1; i<=10; i++)); do
    ffmpeg   -re -i /home/kienpt/Documents/Beam/data/enhanced_1588047114.wav -f s16le -ar 16000 -ac 1 -vsync 1 - > /home/kienpt/audioFiles/virtmic
    #cat ../data/test.raw  > /home/kienpt/audioFiles/virtmic
done