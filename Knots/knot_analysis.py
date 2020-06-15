from array import array
from ROOT import TCanvas, TPad, TTree, TFormula, TF1, TPaveLabel, TH1F, TH2F, TFile
from ROOT import gROOT, gBenchmark, gDirectory,gStyle, TGraph, TLine
from math import sin,cos
from knots import findcenter

gROOT.Reset()

def getCWH(x):
    xx = []
    yy = []
    
    for i in x.values():
        i[0][0].append(x)
        i[1][0].append(x)
        i[0][1].append(y)
        i[1][1].append(y)
    


def Draw_Knot(knot):
    n =len(knot['paths'])
    lines = [TLine(knot['paths'][i][0][0],knot['paths'][i][0][1],knot['paths'][i][1][0],knot['paths'][i][1][1]) for i in range(n)]
    
    c1 = TCanvas() 
    c1.SetFillColor( 29 )

    for i in knot['elevators'].keys():
        if knot['elevators'][i] == 0:
            lines[i].SetLineColor(4)
            lines[i].SetLineWidth(4)
            lines[i].Draw()
        else:
            lines[i].SetLineColor(2)
            lines[i].SetLineWidth(4)

    for i in knot['elevators'].keys():
        if knot['elevators'][i] == 1:
            lines[i].Draw()    
   
    zwom =input("type a number to continue")
    return lines

#print findcenter(h)

#c1 = TCanvas()
#knot = Draw_Knot(h,c1)
#for i in knot:
#    i.Draw()
#zwom = input("sheat")



