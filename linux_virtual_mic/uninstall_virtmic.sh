#!/bin/bash

# Uninstall the virtual microphone.

pactl unload-module module-pipe-source
rm /home/trungkien1102/.config/pulse/client.conf