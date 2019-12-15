#!/bin/sh
# orig: 640x480

ffmpeg -i DropSeq_original.mp4 -filter:v "crop=260:480:190:0,transpose=2" DropSeq.mp4
ffmpeg -i DropSeq.mp4 -r 1 -s 480x260 -f image2 foo-%03d.jpg
ffmpeg -i DropSeq.mp4 -vcodec copy -acodec copy DropSeq.avi
mv foo-001.jpg DropSeq_thumbnail.jpg
rm foo*

