import ROOT as rt
import random

mem=rt.TH1F("mem","mem",100,0,1)
for i in range(10000):
    mem.Fill(random.random())

mem.Draw("E0")
zwom=input("sheat")
