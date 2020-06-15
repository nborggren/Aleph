from MakeBraid import *
import os

print '?!'*21
print 'Welcome to the sBraid production line'
print 'Nathan Borggren'
print 'August 2009'
print '?!'*21

z = raw_input("How many magnetic lines do you want?  ")
mlines = int(z)



print " "
print "The Rules of the Game:"
print "alphabet:"
print "    \'mn\' crosses line n, n+1"
print "    \'Mn\' is the inverse of mn, crossing line n+1 over n"
print "    \'xn'\' makes a copy of line n and switches line n+1 to n XOR n+1" 
print "    \'Xn'\' makes a copy of line n+1 and switches line n to n XOR n+1" 
print " "


eBraid = raw_input("how should I braid them?  ")

l = len(eBraid)/2

allbraids = [eBraid[2*i:2*i+2] for i in range(l)]

x = 2
b = []
ib = []
bb = []
ibb = []

elem = ['braid','vbraid','magbraid','vmagbraid','ibraid','magibraid','vibraid','vmagibraid','bbraid','cbraid','magnot','maginot','magibbraid','dbraid','ident','magident']

gens = {}
gens2 = {}

alph = ['a','A','b','B','c','C','d','D','e','E','f','F','g','G','h','H','i','I']

for i in range(len(elem)):
    gens.update({elem[i]:[]})
    gens2.update({alph[i]:elem[i]})

for i in allbraids:
    print i
    gens[gens2[i[0]]].append([x,int(i[1])])
#     if i[0] == 'm':
#         b.append([x,int(i[1])])
#     elif i[0] == 'M':
#         ib.append([x,int(i[1])])
#     elif i[0] == 'x':
#         bb.append([x,int(i[1])])
#     elif i[0] == 'X':
#         ibb.append([x,int(i[1])])
    x = x+1



    


ids = GetIdsM(x,mlines,gens)

gens.update({'magident':ids})

MakeBraidM("./from_BraidCalc/"+eBraid+".ps",gens)

#os.system("gimp ./from_BraidCalc/"+eBraid+".ps")
