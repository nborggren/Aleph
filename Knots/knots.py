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
    tangles = {0:[x,y]}
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

    elevators = {0:0,1:1}
    players = {0:0,1:1}
    temp = {'points':points,
            'tangles':tangles,
            'paths':paths,
            'elevators':elevators,
            'players':players,
            'corners':corners}
    return temp

#planar translations

def translate(knot,x,y):
    nknot = {'paths':{},'points':{},'corners':{},
            'elevators':{},
             'tangles':{},
            'players':{}}
    for i in ['points','corners','tangles']:
        q=0
        for j in knot[i].keys():
            #print j[0]+x,j[1]+y
            nknot[i].update({j:[knot[i][j][0]+x,knot[i][j][1]+y]})
            q=q+1
    
    q = 0
    for i in knot['paths'].values():
        #print x, i[0][0],i[0][0]+x
        nknot['paths'].update({q:[[i[0][0]+x,i[0][1]+y],[i[1][0]+x,i[1][1]+y]]})
        
        q=q+1
    
    for i in ['elevators','players']:
        for j in knot[i].keys():
            nknot[i].update({j:knot[i][j]})
            
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
    tangles = {}

    temp = {'points':points,
            'paths':paths,
            'elevators':elevators,
            'players':players,
            'tangles':tangles,
            'corners':corners}

    for i in ['elevators','players','tangles']:
        for j in x[i].keys():
            #print x[i][j], j,(1+x[i][j])%2, 'elevate'
            #temp[i].update({j:(1+x[i][j])%2})  #confused
            temp[i].update({j:x[i][j]})  #confused

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
    tangles = {}

    temp = {'points':points,
            'paths':paths,
            'elevators':elevators,
            'players':players,
            'tangles':tangles,
            'corners':corners}

    for i in ['elevators','players','tangles']:
        for j in x[i].keys():
            temp[i].update({j:x[i][j]})


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
    tangles = {}

    temp = {'points':points,
            'paths':paths,
            'elevators':elevators,
            'players':players,
            'tangles':tangles,
            'corners':corners}

    for i in ['elevators','players','tangles']:
        for j in x[i].keys():
            temp[i].update({j:x[i][j]})
    
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
    tangles = {}

    temp = {'points':points,
            'paths':paths,
            'elevators':elevators,
            'players':players,
            'tangles':tangles,
            'corners':corners}

    for i in ['elevators','players','tangles']:
        for j in x[i].keys():
            temp[i].update({j:x[i][j]})

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
    tangles = {}
    
    temp = {'points':points,
            'paths':paths,
            'elevators':elevators,
            'players':players,
            'tangles':tangles,
            'corners':corners}

    for i in ['elevators','players','tangles']:
        for j in x[i].keys():
            temp[i].update({j:x[i][j]})

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
    tangles = {}

    
    temp = {'points':points,
            'paths':paths,
            'elevators':elevators,
            'players':players,
            'tangles':tangles,
            'corners':corners}

    for i in ['elevators','players','tangles']:
        for j in x[i].keys():
            temp[i].update({j:x[i][j]})

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
    tangles = {}

    
    temp = {'points':points,
            'paths':paths,
            'elevators':elevators,
            'players':players,
            'tangles':tangles,
            'corners':corners}

    for i in ['elevators','players','tangles']:
        for j in x[i].keys():
            temp[i].update({j:x[i][j]})

    return temp  

# a certain sense of lack of completeness inspires the last:

def identity(x):
    return x

#end D4 group
#???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
#operations to construct knots

def add(x,y):

    t1x = widthandheight(x)
    t1y = widthandheight(y)
    t2x = findcenter(x)
    t2y = findcenter(y)

    rofx2lofy = (t2y[0]-t1y[0]/2.) - (t2x[0]+t1x[0]/2.)
    #tofx2bofy = (t2y[1]-t1y[1]/2.) - (t2x[1]+t1x[1]/2.)

    if rofx2lofy < 0:
        y = translate(y,-1.25*rofx2lofy,0)

    #if tofx2bofy < 0:
    #    y = translate(y,0,-tofx2bofy)
        
                                          

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

    tangles = {}
    i = 0
    for k in x['tangles'].values():
        points.update({i:k})
        i=i+1
    for k in y['tangles'].values():
        points.update({i:k})
        i=i+1
    
    npoints = len(x['paths'].keys())+len(y['paths'].keys())
    #print npoints
    paths = {}
    elevators = {}
    players = {}

    i = 0
    for k in x['paths'].keys():
        paths.update({i:x['paths'][k]})
        elevators.update({i:x['elevators'][k]})
        players.update({i:x['players'][k]})
        i=i+1
    for k in y['paths'].keys():
        paths.update({i:y['paths'][k]})
        elevators.update({i:y['elevators'][k]})
        players.update({i:y['players'][k]})
        i=i+1

    paths.update({i:[x['corners']['ne'],y['corners']['nw']]})
    elevators.update({i:0})
    players.update({i:0})
    i=i+1
    elevators.update({i:1})
    players.update({i:1})
    paths.update({i:[x['corners']['se'],y['corners']['sw']]})
    corners = {'ne':y['corners']['ne'],
               'nw':x['corners']['nw'],
               'se':y['corners']['se'],
               'sw':x['corners']['sw']}
    

    temp = {'points':points,
            'paths':paths,
            'elevators':elevators,
            'tangles':tangles,
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
    c = ymx(x)
    b = findcenter(c)
    d = findcenter(y)
    #print b
    temp = add(c,translate(y,0,b[1]-d[1]))
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
    for i in x['points'].values():
        #print i
        xs.append(i[0])
        ys.append(i[1])
    xx = min(xs)+(max(xs)-min(xs))/2.
    yy = min(ys)+(max(ys)-min(ys))/2.
    
    return [xx,yy]

def move2origin(x):
    zz = findcenter(x)
    #print 'how far to origin', zz
    return translate(x,-zz[0],-zz[1])

def widthandheight(x):

    xs = []
    ys = []
    for i in x['paths'].values():
        xs.append(i[0][0])
        xs.append(i[1][0])
        ys.append(i[0][1])
        ys.append(i[1][1])
    xx = max(xs)-min(xs)
    yy = max(ys)-min(ys)
    
    return [xx,yy]

def scale(x,aaa):
    knot = move2origin(x)
    b = widthandheight(knot)
    points = {}
    corners = {}
    paths = {}
    nknot = {'paths':{},'points':{},'corners':{},'elevators':{},'players':{},'tangles':{}}
    for i in ['points','corners','tangles']:
        q=0
        for j in knot[i].keys():
            #print j[0]/b[0],j[1]/b[1]
            nknot[i].update({j:[knot[i][j][0]*aaa/b[0],aaa*knot[i][j][1]/b[1]]})
            q=q+1
    q=0    
    for i in knot['paths'].values():
        nknot['paths'].update({q:[[aaa*i[0][0]/b[0],aaa*i[0][1]/b[1]],[aaa*i[1][0]/b[0],aaa*i[1][1]/b[1]]]})
        q=q+1

    for i in ['elevators','players']:
        for j in x[i].keys():
            nknot[i].update({j:x[i][j]})

    nknot = translate(nknot,.5,.5)
            
    return nknot   

def tie(x):
    b = findcenter(x)
    c = widthandheight(x)
    d = [b[0],b[1]+3*c[1]/4.]
    f = [b[0],b[1]-3*c[1]/4.]

    n = len(x['paths'])
    x['paths'].update({n:[x['corners']['ne'],d]})
    x['paths'].update({n+1:[d,x['corners']['nw']]})
    x['paths'].update({n+2:[x['corners']['se'],f]})
    x['paths'].update({n+3:[f,x['corners']['sw']]})
    for i in ['elevators','players']:
        x[i].update({n:0})
        x[i].update({n+1:0})
        x[i].update({n+2:1})
        x[i].update({n+3:1})
    return x
    
#???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????


z = sum([tangle(2*i,0) for i in range(5)])
q = sum([tangle(10+2*i,0) for i in range(5)])
pp = sum([tangle(13+2*i,0) for i in range(3)])
zz = sum([tangle(20+2*i,0) for i in range(7)])

h = product([z,q,pp,zz])



h = tie(h)

h = scale(h,.9)

#f = open('5537.txt','w')
#f.write(str(h))





