
def MakeBraid(name,ids,braids,ibraids):
    f = open("template.ps")
    g = open(name, 'w')
    
    for line in f:
        if line.find("REPLACE ME")== -1:
            g.write(line)
        else:
            for j in ids:
                g.write(str(j[0])+" "+str(j[1])+ "  ident \n")
            for j in braids:
                g.write(str(j[0])+" "+str(j[1])+ "  braid  \n")
            for j in ibraids:
                g.write(str(j[0])+" "+str(j[1])+ "  ibraid  \n")
    g.close()

def MakeBraidM(name,Gens):
    f = open("template.ps")
    g = open(name, 'w')
    
    for line in f:
        if line.find("REPLACE ME")== -1:
            g.write(line)
        else:
            for k in Gens.keys():
                for j in Gens[k]:
                    g.write(str(j[0])+" "+str(j[1])+ "  "+ k + '\n')
    g.close()

def GetIds(n,m,braids,ibraids):
    ids = []
    for i in range(n):
        for j in range(m):
            ids.append([i+1,j+1])

    for j in braids:
        b = [e for i, e in enumerate(braids) if e==j]
        for k in b:
            ids.remove(k)
            ids.remove([k[0],k[1]+1])
            print k,[k[0],k[1]+1] 


    for j in ibraids:
        b = [e for i, e in enumerate(ibraids) if e==j]
        for k in b:
            ids.remove(k)
            ids.remove([k[0],k[1]+1])

    return ids


def GetIdsM(n,m,Gens):
    ids = []
    for i in range(n):
        for j in range(m):
            ids.append([i+1,j+1])

    for elem in Gens.keys():
        for j in Gens[elem]:
            b = [e for i, e in enumerate(Gens[elem]) if e==j]
            for k in b:
                try:
                    if elem in ['magnot','maginot']:
                        ids.remove(k)
                    else:
                        ids.remove(k)
                        ids.remove([k[0],k[1]+1])
                except ValueError:
                    print k,[k[0],k[1]+1]

    return ids

def GetVIdsM(n,m,Gens):
    ids = []
    for i in range(n):
        for j in range(m):
            ids.append([i+1,j+1])

    for elem in Gens.keys():
        for j in Gens[elem]:
            b = [e for i, e in enumerate(Gens[elem]) if e==j]
            for k in b:
                try:
                    ids.remove(k)
                    ids.remove([k[0]+1,k[1]])
                except ValueError:
                    print k,[k[0],k[1]+1]

    return ids

#braids = [[3,4],[2,7]]
#ibraids = [[6,5]]
#dbraids = [[1,9],[3,9],[1,5]]
#bbraids = [[11,3],[10,5]]
#cbraids = [[2,1],[5,1]]

# n = 4
# m = 5

# vmagbraids = []
# vmagibraids = []
# magbraids = []
# magibraids = []
# braids = []
# ibraids = []
# vbraids = []
# vibraids = []

# for i in range(n):
#     braids.append([4*i+2,2*i+2])
#     ibraids.append([4*i+4,2*i+3])
#     braids.append([4*i+4,2*i+1])
#     ibraids.append([4*i+6,2*i+2])
#     vbraids.append([2*i+2,4*i+2])
#     vibraids.append([2*i+3,4*i+4])
#     vbraids.append([2*i+1,4*i+4])
#     vibraids.append([2*i+2,4*i+6])
    

# Gens = {'bbraid':vbraids,
#         'cbraid':vibraids
#         }


# Gens2 = {'magbraid':braids,
#         'magbbraid':ibraids
#         }

# #ids = GetIds(19,19,braids,ibraids)

# print Gens

# ids2 = GetVIdsM(19,19,Gens)
# ids3 = GetIdsM(19,19,Gens2)

# #MakeBraid("test.ps",ids,braids,ibraids)

# Gens.update({'vmagident':ids2})
# for i in Gens2.keys():
#     Gens.update({i:Gens2[i]})
# Gens.update({'magident':ids3})

# MakeBraidM("ring4.ps",Gens)




# #print 'ids', ids
# print 'ids2', ids2

# #zwom = input("sheat")

