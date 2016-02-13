import sys

if len(sys.argv)<2:
    from lorenz import *
    tmp=lorenz()
else:
    print sys.argv[1]
    exec("from "+sys.argv[1]+" import *") 
    exec("tmp = "+sys.argv[1]+"()")

tmp.Init()
tmp.Start()

zwom = input("sheat")
