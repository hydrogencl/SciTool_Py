import os, sys, math, re
import argparse
from pytube import YouTube
from subprocess import run

parser = argparse.ArgumentParser(description='Reading input date')


parser.add_argument("-u","--url",type=str,dest='strURL', 
                    help='url of the youtube, full')

parser.add_argument("-o","--out",type=str,dest='strFMT', default="mp4" ,
                    help='output format of the file')

parser.add_argument("-a","--audio-only", action='store_true', dest='ifAudio', default=False,
                    help='if get only audio files')

#parser.add_argument("-n","--num-o-only", action='store_true', dest='ifAudio', 
#                    help='if get only audio files')

args = parser.parse_args()
print("reading the url ... ")
url  = args.strURL


if args.ifAudio:
    Stream_obj = YouTube(url).streams.filter()
    for obj in Stream_obj:
        #print(obj.abr)
        if obj.abr == '160kbps':
            itag_in = obj.itag
    Stream_obj_out = Stream_obj.get_by_itag(itag_in)
    strFile = Stream_obj_out.download()
    arrFileString = re.split(".", strFile)
    strFileNew = "{0:s}.{1:s}".format(strFile, 'mp3')
    try:
        run(["ffmpeg","-i",strFile, "-ab", "320k", strFileNew])
    except:
        print("Sorry not found! Fail")

