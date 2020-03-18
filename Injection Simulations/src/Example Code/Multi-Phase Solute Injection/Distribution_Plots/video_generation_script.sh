ffmpeg -r 24 -i %04d.png -vcodec libx264 -y -an output_video.mp4 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2"

