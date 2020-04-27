#!/bin/bash

# Uninstall the virtual microphone.

pactl unload-module module-pipe-source
rm /home/kienpt/.config/pulse/client.conf