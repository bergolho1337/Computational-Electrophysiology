#!/bin/bash

# Variables
FILENAME="frames/discordant_alternans_beeler_reuter"
FRAME_RATE="30"
END_FRAME="1240"
OUTPUT_VIDEO_FILENAME="video/discordant_alternans_beeler_reuter"
RESOLUTION="1020x430"

# Execute the converting command using FFMPEG
ffmpeg -r $FRAME_RATE -f image2 -s $RESOLUTION -start_number 1 -i $FILENAME.%04d.jpg -vframes $END_FRAME -vcodec libx264 -crf 25  -pix_fmt yuv420p $OUTPUT_VIDEO_FILENAME.mp4
