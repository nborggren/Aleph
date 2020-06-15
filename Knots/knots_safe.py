#???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
#???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
#
#knots.py
#N.A. Borggren June 2009
#
#designs a simple geometrical realization of a knot topology so as to be embedded in a Cadence layout upon 
#manipulation in accordance with Hypres design rules.
#
#???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????]
#???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????

layers = ['m1', 'm2', 'm3']


#knots formed by products of tangles in conway notation, digits add to # of crossings.  
#Here are the knots up to eight crossings that don't require more sophisticated notation

aknots = [3,22,5,32,42,312,2112,7,52,43,322,313,2212,21112,62,512,44,413,4112,332,3212,3113,31112,2312,2222,22112]

dknots = {}


#the basic tangle unit cell.  Centered at x,y it occupies a 1 by 1 area
def tangle(x,y):
    points = {0:[x+.5,y+.5],
            1:[x-.5,y+.5],
            2:[x-.5,y-.5],
            3:[x+.5,y-.5]}
    corners = {'ne':[x+.5,y+.5],
            'nw':[x-.5,y+.5],
            'sw':[x-.5,y-.5],
            'se':[x+.5,y-.5]}

    paths = {0:[points[1],points[3]],
             1:[points[0],points[2]],
             }

    elevators = {}
    players = {}
    temp = {'points':points,
            'paths':paths,
            'elevators':elevators,
            'players':players,
            'corners':corners}
    return temp

#planar translations

def translate(knot,x,y):
    nknot = knot
    for i in ['points','corners']:
        for j in knot[i].keys():
            nknot[i][j][0] = knot[i][j][0] + x
            nknot[i][j][1] = knot[i][j][1] + y
    for i in knot['paths'].keys():
        for j in range(len(knot['paths'][i])):
            
            nknot['paths'][i][j][0] = knot['paths'][i][j][0] + x
            nknot['paths'][i][j][1] = knot['paths'][i][j][1] + y
    return nknot
            
#???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
#low budget reflections   
#multiplication requires first a reflection of the first knot in the line y = -x (ymx)

def ymx(x):
    points = {}
    #print 'in spin', x['points']
    for i in x['points'].keys():
        points.update({i:[-x['points'][i][1],-x['points'][i][0]]})
    paths = {}
    for i in x['paths'].keys():
        paths.update({i:[[-x['paths'][i][0][1],-x['paths'][i][0][0]],[-x['paths'][i][1][1],-x['paths'][i][1][0]]]})
    corners = {'ne':[-x['corners']['sw'][1],-x['corners']['sw'][0]],
               'nw':[-x['corners']['nw'][1],-x['corners']['nw'][0]],
               'se':[-x['corners']['se'][1],-x['corners']['se'][0]],
               'sw':[-x['corners']['ne'][1],-x['corners']['ne'][0]]}
    
    elevators = {}
    players = {}
    temp = {'points':points,
            'paths':paths,
            'elevators':elevators,
            'players':players,
            'corners':corners}
    #print 'in spin', temp
    return temp  

def ypx(x):
    points = {}
    #print 'in spin', x['points']
    for i in x['points'].keys():
        points.update({i:[x['points'][i][1],x['points'][i][0]]})
    paths = {}
    for i in x['paths'].keys():
        paths.update({i:[[x['paths'][i][0][1],x['paths'][i][0][0]],[x['paths'][i][1][1],x['paths'][i][1][0]]]})
    corners = {'ne':[x['corners']['sw'][1],x['corners']['sw'][0]],
               'nw':[x['corners']['nw'][1],x['corners']['nw'][0]],
               'se':[x['corners']['se'][1],x['corners']['se'][0]],
               'sw':[x['corners']['ne'][1],x['corners']['ne'][0]]}
    
    elevators = {}
    players = {}
    temp = {'points':points,
            'paths':paths,
            'elevators':elevators,
            'players':players,
            'corners':corners}
    #print 'in spin', temp
    return temp  

#reflection about y is zero x
def yi0x(x):
    points = {}
    #print 'in spin', x['points']
    for i in x['points'].keys():
        points.update({i:[x['points'][i][0],-x['points'][i][1]]})
    paths = {}
    for i in x['paths'].keys():
        paths.update({i:[[x['paths'][i][0][0],-x['paths'][i][0][1]],[x['paths'][i][1][0],-x['paths'][i][1][1]]]})
    corners = {'ne':[x['corners']['sw'][0],x['corners']['sw'][1]],
               'nw':[x['corners']['nw'][0],x['corners']['nw'][1]],
               'se':[x['corners']['se'][0],x['corners']['se'][1]],
               'sw':[x['corners']['ne'][0],x['corners']['ne'][1]]}
    
    elevators = {}
    players = {}
    temp = {'points':points,
            'paths':paths,
            'elevators':elevators,
            'players':players,
            'corners':corners}
    #print 'in spin', temp
    return temp  


#reflection about y is infinity (omega) x
def yiwx(x):
    points = {}
    #print 'in spin', x['points']
    for i in x['points'].keys():
        points.update({i:[-x['points'][i][0],x['points'][i][1]]})
    paths = {}
    for i in x['paths'].keys():
        paths.update({i:[[-x['paths'][i][0][0],x['paths'][i][0][1]],[-x['paths'][i][1][0],x['paths'][i][1][1]]]})
    corners = {'ne':[-x['corners']['sw'][0],x['corners']['sw'][1]],
               'nw':[-x['corners']['nw'][0],x['corners']['nw'][1]],
               'se':[-x['corners']['se'][0],x['corners']['se'][1]],
               'sw':[-x['corners']['ne'][0],x['corners']['ne'][1]]}
    
    elevators = {}
    players = {}
    temp = {'points':points,
            'paths':paths,
            'elevators':elevators,
            'players':players,
            'corners':corners}
    #print 'in spin', temp
    return temp  

def piby2(x):
    points = {}
    #print 'in spin', x['points']
    for i in x['points'].keys():
        points.update({i:[-x['points'][i][1],x['points'][i][0]]})
    paths = {}
    for i in x['paths'].keys():
        paths.update({i:[[-x['paths'][i][0][1],x['paths'][i][0][0]],[-x['paths'][i][1][1],x['paths'][i][1][0]]]})
    corners = {'ne':[-x['corners']['sw'][1],x['corners']['sw'][0]],
               'nw':[-x['corners']['nw'][1],x['corners']['nw'][0]],
               'se':[-x['corners']['se'][1],x['corners']['se'][0]],
               'sw':[-x['corners']['ne'][1],x['corners']['ne'][0]]}
    
    elevators = {}
    players = {}
    temp = {'points':points,
            'paths':paths,
            'elevators':elevators,
            'players':players,
            'corners':corners}
    #print 'in spin', temp
    return temp  



#rotation by pi
def pi(x):
    points = {}
    #print 'in spin', x['points']
    for i in x['points'].keys():
        points.update({i:[-x['points'][i][0],-x['points'][i][1]]})
    paths = {}
    for i in x['paths'].keys():
        paths.update({i:[[-x['paths'][i][0][0],-x['paths'][i][0][1]],[-x['paths'][i][1][0],-x['paths'][i][1][1]]]})
    corners = {'ne':[-x['corners']['sw'][0],-x['corners']['sw'][1]],
               'nw':[-x['corners']['nw'][0],-x['corners']['nw'][1]],
               'se':[-x['corners']['se'][0],-x['corners']['se'][1]],
               'sw':[-x['corners']['ne'][0],-x['corners']['ne'][1]]}
    
    elevators = {}
    players = {}
    temp = {'points':points,
            'paths':paths,
            'elevators':elevators,
            'players':players,
            'corners':corners}
    #print 'in spin', temp
    return temp  

def pi3by2(x):
    points = {}
    #print 'in spin', x['points']
    for i in x['points'].keys():
        points.update({i:[-x['points'][i][1],x['points'][i][0]]})
    paths = {}
    for i in x['paths'].keys():
        paths.update({i:[[-x['paths'][i][0][1],x['paths'][i][0][0]],[-x['paths'][i][1][1],x['paths'][i][1][0]]]})
    corners = {'ne':[-x['corners']['sw'][1],x['corners']['sw'][0]],
               'nw':[-x['corners']['nw'][1],x['corners']['nw'][0]],
               'se':[-x['corners']['se'][1],x['corners']['se'][0]],
               'sw':[-x['corners']['ne'][1],x['corners']['ne'][0]]}
    
    elevators = {}
    players = {}
    temp = {'points':points,
            'paths':paths,
            'elevators':elevators,
            'players':players,
            'corners':corners}
    #print 'in spin', temp
    return temp  

# a certain sense of lack of completeness inspires the last:

def identity(x):
    return x

#end D4 group
#???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
#operations to construct knots

def add(x,y):
    npoints = len(x['points'].keys())+len(y['points'].keys())
    #print npoints
    points = {}
    i = 0
    for k in x['points'].values():
        points.update({i:k})
        i=i+1
    for k in y['points'].values():
        points.update({i:k})
        i=i+1
    
    npoints = len(x['paths'].keys())+len(y['paths'].keys())
    #print npoints
    paths = {}
    i = 0
    for k in x['paths'].values():
        paths.update({i:k})
        i=i+1
    for k in y['paths'].values():
        paths.update({i:k})
        i=i+1    
    paths.update({i:[x['corners']['ne'],y['corners']['nw']]})
    paths.update({i:[x['corners']['se'],y['corners']['sw']]})
    corners = {'ne':y['corners']['ne'],
               'nw':x['corners']['nw'],
               'se':y['corners']['se'],
               'sw':x['corners']['sw']}
    
    elevators = {}
    players = {}
    temp = {'points':points,
            'paths':paths,
            'elevators':elevators,
            'players':players,
            'corners':corners}
    return temp

#sum is to add an array of tangles left to right
def sum(tangs):
    temp = add(tangs.pop(0),tangs.pop(0))
    while len(tangs)>0:

        temp = add(temp, tangs.pop(0))

    return temp

def multiply(x,y):
    temp = add(ymx(x),y)
    return temp

def product(tangs):
    temp = multiply(tangs.pop(0),tangs.pop(0))
    
    while len(tangs)>0:
        
        temp = multiply(temp, tangs.pop(0))
    
    return temp            

#???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
#some circuit maintenance functions

def findcenter(x):
    xs = []
    ys = []
    for i in x['points']:
        xs.append(i[0])
        ys.append(i[1])
    xx = (max(xs)-min(xs))/2.
    yy = (max(ys)-min(ys))/2.
    
    return [xx,yy]

def move2origin(x):
    zz = findcenter(x)
    return translate(x,-zz[0],-zz[1])
    
#???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????



a = tangle(0,0)
b = tangle(0,1)
c = tangle(0,2)

d = add(a,b)
e = add(d,c)
f = sum([a,b,c])

g = translate(d,2,1)
print 'd', d
print 'e', e
print 'f', f
print 'g', g

#print 'spin', spin(d)
#print 'multiply', multiply(a,b)
print len(aknots), 'aknots'
print 'product', product([a,b,c])

z = multiply(f,g)
print '32?'

for i in z['points'].values():
    print i

for i in z['paths'].values():
    print i
