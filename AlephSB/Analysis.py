from dolfin import *
from os import system as MySys
from ROOT import TFile, gStyle, gDirectory, TCanvas, TGraph2D, TGraph
from array import array as MyArray
from os import environ
##@package Analysis 
#A class for investigating the dolfin solutions.

## the data count as a len d string with zeroes on the left followed by the stringed int for use in naming the indexed output files. 
#
#e.g 
#
# >>> from Analysis import *
#
# >>> count=cname(0)
#
# >>> for i in range(3):count.next()
#
# '0000000'
#
# '0000001'
#
# '0000002'
#
# >>> a=cname(0,D=14)
#
# >>> a.next()
#
# '00000000000000'
def cname(i,D=6):
    while 0<1:
        ni = len(str(i))
        si = (D-ni)*'0'+str(i)
        i=i+1
        yield si

def PostAnalysis(name,rng,integrals = ["u1-u0"]):
    cnt = cname(1)
    mesh = Mesh(name+"_mesh.xml")
    Q = FunctionSpace(mesh,"Lagrange",2)
    u = Function(Q)
    V = VectorFunctionSpace(mesh,"Lagrange",2)
    v = Function(V)

    for i in range(rng):
        tmp = cnt.next()
        if i == 0:
            u1 = Function(Q,name+"_"+tmp+".xml.gz")
        else:
            u0 = u1
            u1 = Function(Q,name+"_"+tmp+".xml.gz")
        u.interpolate(u1)
        p=plot(u)


## finds the difference between U0,U1, and returns in du
# (rather project(U0-U1)) homework: show this functions utility; rates of convergence
def Du(U0,U1):
    return project(U0-U1)

## the entropy density evaluates the entropy, \f$-u ln(u)\f$, of a function u at the mesh coordinates and puts into a ROOT TGraph2D object.  homework: extend for D>2, analytically continue into complex domain 
def entropy_density(mesh,u,name="test",D=2):
    x,y,z = MyArray('f'),MyArray('f'),MyArray('f')
    for i in mesh.coordinates()[:]:
        tmp = sqrt(pow(u(i),2))
        x.append(i[0])
        y.append(i[1])
        z.append(-tmp*ln(tmp))
    return TGraph2D(len(x),x,y,z)

## the entropy density array zips through X,Y and evaluates the entropy, \f$-u ln(u)\f$, of a function u at those coordinates and puts into a ROOT TGraph object.  homework: extend for D>2, analytically continue into complex domain, (postprocessing of TGraphs, TGraph2Ds)
def entropy_density_array(X,Y,u,name="test",D=2,t=0):

    n=len(X)
    if t==0:
        t=[float(i)/n for i in range(n)]

    jj=0
    x,z = MyArray('f'),MyArray('f')
    for i,j in zip(X,Y):
        tmp = sqrt(pow(u([i,j]),2))
        x.append(t[jj])
        z.append(-tmp*ln(tmp))
        jj+=1
    return TGraph(n,x,z)

## like entropy, but for the potential, \f$- ln(u)\f$, outputs ROOT TGraph2D 
def energy_density(mesh,u):

    x,y,z = MyArray('f'),MyArray('f'),MyArray('f')
    for i in mesh.coordinates()[:]:
        tmp = sqrt(pow(u(i),2))
        x.append(i[0])
        y.append(i[1])
        z.append(-ln(tmp))
    tmp = TGraph2D(len(x),x,y,z)
    return tmp

## like entropy, but for the potential, \f$- ln(u)\f$, outputs ROOT TGraph
def energy_density_array(X,Y,u,name="test",D=2,t=0):
    n=len(X)
   
    if t==0:
        t=[float(i)/n for i in range(n)]

    jj=0
    x,z = MyArray('f'),MyArray('f')
    for i,j in zip(X,Y):
        tmp = sqrt(pow(u([i,j]),2))
        x.append(t[jj])
        z.append(-ln(tmp))
        jj+=1

    return TGraph(n,x,z)

## some day quality will be checked
def QualCheck(mesh,fname,D,vel):
    pass

## quick load for evaluating, returns dolfin mesh, and function space
def Load_2D(n,d=2,x=[0,1],y=[0,1],elements = "Lagrange"):
    mesh=Rectangle(x[0],y[0],x[1],y[1],n,n)
    Q=FunctionSpace(mesh,elements,d)
    return mesh,Q

## quick load for evaluating say, the lorenz poincare section
def Load_3D(name,mins=[-30,-30,26],maxs=[30,30,28],nx=30,ny=30,nz=4):
    mesh = Box(mins[0],mins[1],mins[2],maxs[0],maxs[1],maxs[2],nx,ny,nz)
    Q=FunctionSpace(mesh,"Lagrange",2)
    u = Function(Q,name)
    return mesh,Q,u

## quick loads the file named u as well for evaluating
def Load_Func(u,n,d=2,x=[0,1000],y=[0,1000],elements = "Lagrange"):
    mesh,Q = Load_2D(n,d=d,x=x,y=y,elements = elements)
    return mesh,Q,Function(Q,u)
    
## the q derivative entails a higher density of grid points as one approaches zero for qx,qy < 1, in accordance with \f$x_out=q_x^{x_in*n_x}\f$.  It approaches the ordinary derivative as qx,qy approach 1.
def qUnitSquare(nx,ny,qx,qy):
    mesh = UnitSquare(nx,ny)
    for x in mesh.coordinates():
        x[0]=qx**(x[0]*nx)
        x[1]=qy**(x[1]*ny)
    return mesh

## returns a TGraph2D of the constant z cross section (z keyword) of a 3D function u given the mesh keyword.
def Lorenz_Psect(u,mesh = Rectangle(-30,-30,30,30,40,40),z=27):

    x,y,z = MyArray('f'),MyArray('f'),MyArray('f')
    for i in mesh.coordinates()[:]:
            x.append(i[0])
            y.append(i[1])
            z.append(u([i[0],i[1],z]))
    return TGraph2D(len(x),x,y,z)

## the q derivative entails a higher density of grid points as one approaches zero for qx,qy < 1, in accordance with \f$x_out=q_x^{x_in*n_x}*(x_max-x_min)\f$.  It approaches the ordinary derivative as qx,qy approach 1.
def qRectangle(nx,ny,qx,qy,mins=[0,0],maxs=[1,1]):
    mesh = UnitSquare(nx,ny)
    for z in mesh.coordinates():
        z[0]=qx**(z[0]*nx)*(maxs[0]-mins[0])
        z[1]=qy**(z[1]*ny)*(maxs[1]-mins[1])
    return mesh

## the qs and hs will take over the world
def qhUnitSquare(hx,hy,qx,qy):
    nx,ny=int(1/hx),int(1/hy)
    mesh = UnitSquare(nx,ny)

    for x in mesh.coordinates():
        ix = x[0]*nx
        iy = x[1]*ny
        x[0] = qh(1,ix,qx,hx)
        x[1] = qh(1,iy,qy,hy)
    return mesh

## shen?
def qh(x,i,q,h):
    if i<=0:
        return x
    else:
        return qh(q*x+h,i-1,q,h) 

## shen?
def qpUnitSquare(q2x,q2y,q1x,q1y,q0x,q0y):
    nx,ny=int(1/q0x),int(1/q0y)
    mesh = UnitSquare(nx,ny)

    for x in mesh.coordinates():
        ix = x[0]*nx
        iy = x[1]*ny
        x[0] = qp(1,ix,q2x,q1x,q0x)
        x[1] = qp(1,iy,q2y,q1y,q0y)
    return mesh

##you cant believe everything doxygen outputs
def qp(x,i,q2,q1,q0):
    if i<=0:
        return x
    else:
        return qp(q2*x*x+q1*x+q0,i-1,q2,q1,q0) 

##Makes movies homework:get rid of token use mesh and cname
def MakeMovies(tok,rng,path="./",P=1,D=0,S=0):
    mesh = Mesh(path+tok+".xml")
    Q = FunctionSpace(mesh,"CG",1)
    tmp,du = Function(Q),Function(Q)
    
    for i in range(rng[0],rng[1]):
        z = Function(Q,path+tok+"_"+(4-len(str(i+1)))*'0'+str(i+1)+".xml")
        if P==1:
            tmp.assign(z)
            v = plot(tmp)
            v.write_png(tok+"_"+(4-len(str(i+1)))*'0'+str(i+1)+"_P.png")
        if S==1:
            du = LR(z,du)
            tmp.assign(du)
            v = plot(tmp)
            v.write_png(tok+"_"+(4-len(str(i+1)))*'0'+str(i+1)+"_S.png")
        if D==1 and i>0:
            w = Function(Q,path+tok+"_"+(4-len(str(i)))*'0'+str(i)+".xml")
            du=Du(w,z,du)
            tmp.assign(du)
            v = plot(du)
            v.write_png(tok+"_"+(4-len(str(i+1)))*'0'+str(i+1)+"_D.png")
      
           
    if P==1:
        system("ffmpeg -qscale 5 -r 5 -b 9600 -i "+tok+"_%06d_P.png "+tok+"_P.mp4")
    if S==1:
        system("ffmpeg -qscale 5 -r 5 -b 9600 -i "+tok+"_%06d_S.png "+tok+"_S.mp4")
    if D==1:
        system("ffmpeg -qscale 5 -r 5 -b 9600 -i "+tok+"_%06d_D.png "+tok+"_D.mp4")
  
##A generator equivalent to newtons method
def NewtonsMethod(y,iyprime,x):
    
    try:
        while 1>0:
            x = x-iyprime(x)*y(x)
            yield x
    except IndexError:
        while 1>0:
            tmp = [x[i] for i in len(x)]
            for i in range(len(x)):
                for j in range(len(x)):
                    tmp[i]=tmp[i]-iyprime[i,j](x)*y[j](x)
                x[i]=tmp[i]
            yield x

##An order parameter that assigns to each coordinate in X the value it has at the 
def NewtonsOrder2D(y,iyprime,X,nmax=100):

    x,y,z1,z2 = MyArray('d'),MyArray('d'),MyArray('d'),MyArray('d')

    for i in X:
        tmp = NewtonsMethod(y,iyprime,i)
        for j in range(nmax):
            out = tmp.next()
        x.append(i[0])
        y.append(i[1])
        z1.append(out[0])
        z2.append(out[1])

    return TGraph2D(len(x),x,y,z2),TGraph2D(len(x),x,y,z2)   

##same complaints as MakeMovies, evaluate function instead of vector C-  
def xml2oct(tok,rng,n):
    mesh = Mesh(tok+".xml")
    Q=FunctionSpace(mesh,"CG",1)
    for i in range(rng):
        if i%n==1:
            w = Function(Q,tok+"_"+str(i)+".xml")
            f=open(tok+"_"+str(i)+".dat",'w')
            z=mesh.coordinates()[:]
            for j in range(len(z)):
                f.write(str(z[j,0])+" "*3+str(z[j,1])+" "*3+str(w.vector()[j])+" \n")

def ShrinkImage(figure,ftype=".eps"):
    MySys("gs -r300 -dEPSCrop -dTextAlphaBits=4 -sDEVICE=png16m -sOutputFile="+figure+".png -dBATCH -dNOPAUSE "+figure+ftype) 
    MySys("convert "+figure+".png eps3:"+figure+"_c.eps")  

def ShrinkImages(filelist,ftype=".eps"):
    f=open(filelist)
    for line in f:
        ShrinkImage(line.split(ftype)[0],ftype=ftype)

##generator? welp let us see what it does:
#
#>>> a=RecRelations(1,2,3,4,5)
#
# >>> a.next()
#
# 1
#
# >>> a.next()
#
# 1
#
# >>> a.next()
#
# 26
#
# >>> a.next()
#
# 201
#
# >>> a.next()
#
# 1901
#
# >>> a.next()
#
# 17126
#
# >>> a.next()
#
# 156001
#
#I dont know anything about that sequence but any number in the last spot after 0,1,0,1 is golden. 
#
# >>> b=RecRelations(0,1,0,1,29301299400382)
#
# >>> b.next()
#
# 1
#
# >>> b.next()
#
# 1
#
# >>> b.next()
#
# 2
#
# >>> b.next()
#
# 3
#
# >>> b.next()
#
# 5
#
# >>> b.next()
#
# 8
#
# >>> b.next()
#
# 13
#
# >>> b.next()
#
# 21
#
# >>> b.next()
#
# 34
def RecRelations(a,b,c,d,x):

    fn2 = 1 # "f_{n-2}"
    fn1 = 1 # "f_{n-1}"
    while True:
        (fn1,fn2,oldfn2) = ((a*x+b)*fn1+(c*x+d)*fn2,fn1,fn2)
        yield oldfn2
    
## Make plots from root outputs. where are these root outputs?
def Prettify_Root(MyFile, sHistos): #self,tt,xx,uu=None,store=0):

    f = TFile(MyFile)
    gStyle.SetPalette(1)
    hist = [gDirectory.Get(i) for i in sHistos]
    c2=TCanvas()
    H = len(hist)
    if H%2==0:
        c2.Divide(H/2,2)
    else:
        c2.Divide((H+1)/2,2)
#    c2 = TCanvas()
    
    gColors=[2,4,6,7,8,9,11,12,13,14,15,16,17,18]
    gShapes=[27,26,23,22,28,15,11]

    for i in range(H):
        c2.cd(i+1)
        # hist[i].SetLineColor(4)
        # hist[i].SetLineWidth(1)
        #     grsX[i].SetMarkerColor(4)
        #     grsX[i].SetMarkerStyle(7)
        #     grsX[i].SetTitle("x_{"+str(i)+"}")
        #     grsX[i].GetYaxis().SetTitle("x_{"+str(i)+"}")
        #     grsX[i].GetXaxis().SetTitle("time")
        hist[i].Draw()
    zwom = input("sheat")
    if zwom == 1:
        return hist
    else:
        return

## homework given a set of u variables, predict the most likely x variable.
def histUtoX(MyHists):
    pass

## let us forgive
def T(ti,q,nmax=100,s=False):
    if s == True:
        dq = [q[i]*ti**i for i in range(len(q))]
        return expand(sum(dq))
    else:
        dq = [q[i]*pow(ti,i) for i in range(len(q))]
    
    return sum(dq)

##This gives a generator for a trajectory with initial condition x and time step dt.  Here is an example use for the lorenz attractor. 
#
#>>> from Dynamics import *
#
# >>> from Analysis import *
#
# >>> dyn=Dynamics()
#
# >>> dyn.Init(imagine=1)
#
# >>> x0=dyn.rzero()
#
# >>> xt=qIntegrate(x0,.001,dyn.Flow)
#
# >>> x0
#
# [(-20.698756946446188+2.5449184142103549j), (-27.741584663422792+2.1660923516583601j), (2.6240628528889474+2.1327704439181305j)]
#
# >>> xt.next()
#
# [(-20.769185223615956+2.541130153584835j), (-28.233665707279705+2.2726516460714548j), (3.1857691418927421+2.0116475673078296j)]
#
# >>> xt.next()
#
# [(-20.843830028452594+2.5384453685097013j), (-28.715691540154349+2.3752154656269879j), (3.7578888865040811+1.8873366315339923j)]
#
# >>> xt.next()
#
# [(-20.922548643569613+2.5368130694808739j), (-29.18748339126363+2.4737168487944792j), (4.3403834883879213+1.7599019322097629j)]
#
# >>> xt.next()
#
# [(-21.005197991046554+2.5361821072740098j), (-29.648850843002137+2.5680847901158188j), (4.9332103162237635+1.6294092099590032j)]
#
# >>> xt.next()
#
# [(-21.09163451956611+2.5365011341024277j), (-30.09959199800096+2.6582443476575728j), (5.5363219401834467+1.49592610424901j)]
#
# >>> xt.next()
#
# [(-21.181714094350458+2.5377185662379791j), (-30.539493675345845+2.7441167748438504j), (6.1496653688169625+1.3595226004875669j)]
#
# >>> xt.next()
# [(-21.27529189016041+2.5397825483240379j), (-30.968331636949625+2.8256196769694077j), (6.7731812886967369+1.2202714698165806j)]
#
# >>> xt.next()
# [(-21.372222287628304+2.5426409196104918j), (-31.385870845111913+2.9026671927174275j), (7.4068033073105033+1.0782487009233956j)]
#
def qIntegrate(x,dt,f,data=False):
    xi=x
    while True:
        tmp = [j+f(xi)[i]*dt for i,j in enumerate(xi)]  
        #data['traj'].append(xi)
        xi = [i for i in tmp]
        yield xi

#compress (gz) each file in a list
def comList(fname):
    f = open(fname)
    for line in f:
        MySys("gzip "+line.replace('\n',''))
    return

#the probability current \f$J_i=-D_{ij}\partial_ju+F_iu\f$
def J(V,D,velocity,u):
    return project(-D*grad(u)+velocity*u,V)

# \f$J_i^2f$
def Jsquared(V,D,velocity,u):
    z = J(V,D,velocity,u)
    return dot(z,z)

# \f$\partial_iJ_i\f$
def DivJ(U,V,D,velocity,u):
    z = J(V,D,velocity,u)
    return project(div(z),U)

def CurlJ(U,V,D,velocity,u):
    #diagonal diffusion for starters
    z = J(V,D,velocity,u)
    return project(curl(z),U)

def CurlJ_2D(J, mesh, t=100, eps = 0.001,xdom=[0,1],ydom=[0,1]):
 
    x,y,z=MyArray('f'),MyArray('f'),MyArray('f') 

    for i in mesh.coordinates()[:]:
        x.append(i[0])
        y.append(i[1])
        tmp = 0
        for j in range(t):
            xt = i[0]+eps*sin(2*pi*j/t)
            yt = i[1]+eps*cos(2*pi*j/t)

            #if xt < xdom[0] or xt > xdom[1] or yt < ydom[0] or yt > ydom[1]:
            #    jtmp = [0,0]
            #else:
            jtmp = J([xt,yt])
            #print jtmp

            tmp = tmp+jtmp[0]*eps*cos(2*pi*j/t)-jtmp[1]*eps*sin(2*pi*j/t)
        print tmp
        z.append(tmp/(2*pi*eps))

    return TGraph2D(len(mesh.coordinates()[:]),x,y,z)   

def entropy_prod(Q,head,n,dt):
    t = MyArray('f')
    ent = MyArray('f')
    v = Function(Q)
    for i in range(n):
        ss = (6-len(str(i)))*"0"+str(i)
        g = Function(Q,head+ss+".xml.gz")
        t.append(i*dt)
        v,S = entropy(g,v)
        ent.append(S)
        print S, i
    return v, TGraph(n,t,ent)

def Func2Graph(mesh,func,Z=0,polar=0):

    try:
        MyRange = mesh.coordinates()[:]

    except AttributeError:
        MyRange=mesh
        pass
    
    n,d=len(MyRange),len(MyRange[0])
    
    print n,d

    if d==2:
        try:
            tmp = len(func(MyRange[0]))
        except TypeError:
            tmp=1

        if tmp==1:
            x,y,z=MyArray('f'),MyArray('f'),MyArray('f') 
            for i in MyRange:
                x.append(i[0])
                y.append(i[1])

                if Z==1:
                    z.append(func(i[0],i[1]))  # who cares about units?
                else:
                    z.append(func(i))

            return TGraph2D(n,x,y,z)

        elif tmp==2 and polar==0:
            x,y,zx,zy=MyArray('f'),MyArray('f'),MyArray('f'),MyArray('f') 
            for i in MyRange:
                x.append(i[0])
                y.append(i[1])

                if Z==1:
                    tmp2 = func(i[0],i[1])
                    zx.append(tmp2[0])
                    zy.append(tmp2[1])
                else:
                    tmp2 = func(i)
                    zx.append(tmp2[0])
                    zy.append(tmp2[1])

            return TGraph2D(n,x,y,zx),TGraph2D(n,x,y,zy)

        elif tmp==2 and polar==1:
            x,y,zx,zy=MyArray('f'),MyArray('f'),MyArray('f'),MyArray('f') 
            for i in MyRange:
                x.append(i[0])
                y.append(i[1])

                if Z==1:
                    tmp2 = func(i[0],i[1])
                    zx.append(pow(pow(tmp2[0],2)+pow(tmp2[1],2),.5))
                    zy.append(atan(tmp2[1]/tmp2[0]))
                else:
                    tmp2 = func(i)
                    zx.append(pow(pow(tmp2[0],2)+pow(tmp2[1],2),.5))
                    zy.append(atan(tmp2[1]/tmp2[0]))

            return TGraph2D(n,x,y,zx),TGraph2D(n,x,y,zy)

    else:
        x,y=MyArray('f'),MyArray('f') 
        for i,e in enumerate(MyRange):
            x.append(i)
            y.append(func(e))
        return TGraph(n,x,y)

def OrderParam(mesh,velocity):

    try:
        MyRange = mesh.coordinates()[:]

    except AttributeError:
        MyRange=mesh
        pass
    
    n,d=len(MyRange),len(MyRange[0])
  
    if d==2:
        x,y,z=MyArray('f'),MyArray('f'),MyArray('f') 
        for i in MyRange:
            x.append(i[0])
            y.append(i[1])
            tmp = [j(i) for j in velocity]
            
            if tmp[0]>=0 and tmp[1]>=0:
                z.append(1)
            elif tmp[0]>=0 and tmp[1]<0:
                z.append(0.33)
            elif tmp[0]<0 and tmp[1]>=0:
                z.append(-0.33)
            else:
                z.append(-1)

        return TGraph2D(n,x,y,z)

    else:
        x,y=MyArray('f'),MyArray('f') 
        for i,e in enumerate(MyRange):
            x.append(i)
            y.append(func(e))
        return TGraph(n,x,y)

#try newtons method dork.
def FindZeroes(mesh,func,Z=0,tol=0.001):

    try:
        MyRange = mesh.coordinates()[:]

    except AttributeError:
        MyRange=mesh
        pass
    
    x,y=MyArray('f'),MyArray('f') 

    for i in MyRange:

        if Z==1:
            tmp=func(i[0],i[1])
            if tmp < tol:
                x.append(i[0])
                y.append(i[1])

        else:
            tmp=func(i)
            print func(i),tmp
            if abs(tmp) < tol:
                x.append(i[0])
                y.append(i[1])

    return TGraph(len(x),x,y)


def Entropy2Graph(D,velocity,mesh,Q,FuncList,dt=1,T=310,n=1085):

    n = len(FuncList)
    t,s=MyArray('f'),MyArray('f') 

    for i,j in enumerate(FuncList):
        #j = (6-len(str(i*dn+1)))*"0"+str(i+1)
        u = Function(Q,j)
        v = Function(Q)
        #J = velocity*u-D*grad(u)
        ss = project(inner(T*grad(LR(u,v))-velocity,velocity*u-D*grad(u)))
        t.append(i*dt)
        s.append(assemble(ss*dx,mesh=mesh))
        gr = TGraph(n,t,s)
        gr.SetMarkerColor(4)
        gr.SetMarkerStyle(7)
        gr.SetTitle("entropy production")

    return gr           

def NegEntropyReg():
    return #z = [i for i in mesh.coordinates()[:] if ss(i)<0]

def FunctionArray(Q,head,tail,n):
    return [Function(Q,head+"0"*(6-len(str(i)))+str(i)+tail) for i in range(n)]

def Func2Table(Q,head,n,dn):
    
    g = Function(Q)
    for i in range(n):
        j = (6-len(str(i*dn+1)))*"0"+str(i*dn+1)
        f=Function(Q,head+j+".xml.gz")
        g.interpolate(f)
        p=plot(g,interactive=True)
        p.write_png(head+j+".png")
    
def FuncMax(mesh,func):
    mx = max(func.vector()[:])
    tmp = [i for i, e in enumerate(func.vector()[:]) if e==mx]
    print tmp
    return mesh.coordinates()[tmp[0]]

def CutCroRegion(mesh,func):
    nn=len(mesh.coordinates()[:])
    for i in range(nn):
        if mesh.coordinates()[i][1]>mesh.coordinates()[i][0]:
            func.vector()[i]=0
    return func

def CutCiRegion(mesh,func):
    nn=len(mesh.coordinates()[:])
    for i in range(nn):
        if mesh.coordinates()[i][1]<mesh.coordinates()[i][0]:
            func.vector()[i]=0
    return func

def MyEval(Q,funcs,points):
    
    tmp = []
    for i in funcs:
        tmp.append([])
        f=Function(Q,i)
        for j in points:
            z = f(j)
            tmp[-1].append(z)
            print z
    return tmp

def fpt(mesh,Q,velocity,D):
    return

##partitions a given mesh into n sectors in accordance with a particular output distribution
#4
def Pools(n,fmesh,fname):
    mesh = Mesh(fmesh)
    Q = FunctionSpace(mesh,"CG",1)
    f = Function(Q,fname)
    fmin = min(f.vector()[:])
    fmax = max(f.vector()[:])
    df = (fmax-fmin)/n
    pools = {}
    [pools.update({i:[]}) for i in range(n+1)]
    cds = [fmin+df*i for i in range(n+1)]
    nn = len(f.vector()[:])
    print nn, "thelength"

    for j in range(nn):
        isin = [ i for i, e in enumerate(cds) if f.vector()[j]>=e]
        pools[max(isin)].append(j)
    return pools

##generates a langevinesque trajectory from reservoir exchange within pools. Equipartition states all states of particular probability equally probable, so we can switch to any state on a given probability contour.
#2
def resExch(nn,steps,pools):
    
    x0 = choice(range(nn))
    traj = []
    
    for i in range(steps):
        where = [i for i, e in enumerate(pools.keys()) if x0 in pools[e]]
        print where
        x0 = choice(pools[where[0]])
        traj.append(x0)
                    
    return traj

def Traj2Bin(X,rule):
    sym = []
    for x in X:
        if eval(rule)==True:
            sym.append(1)
        else:
            sym.append(0)
    return sym
    
def DirtyIntegrator(velocity,x,tf,dt): 
    X=[x[0]]
    Y=[x[1]]
    T=[0]
    t=0
    while t<tf:
        tmp=[x[0]+velocity[0](x)*dt,x[1]+velocity[1](x)*dt]
        #print tmp
        x=tmp
        X.append(x[0])
        Y.append(x[1])
        t=t+dt
        T.append(t/60.)
        
    return X,Y,T

def Arrays2Graph(X,Y):
    x=MyArray('f')
    y=MyArray('f')
    for i,j in zip(X,Y):
        x.append(i)
        y.append(j)
    return TGraph(len(x),x,y)
    
def PlotViewer(FileList,names,option="APL"):
    
    for i in FileList:
        f=TFile(i)
        for j in names:
            tmp=gDirectory.Get(j)
            tmp.Draw(option)
            zwom=input("sheat")
    return
        

def LamStates():
    tmp = ['1','2','3']
    return [i+j+k for i in tmp for j in tmp for k in tmp]

#seee TrajTable ~/Analysis/Diff/lambda                  
def TrajTable(mesh,Q,head,nx,dat=10,image=40,zeroes=[],entropy=1,energy=1):
    return

def Publish(name,title,media,columns=5):
    
    f=open(environ["ALEPHPATH"]+'lib/template.tex')
    g=open(name+'.tex','w')

    img_beg = ["\\begin{figure}[h]"
               "\\begin{center}"]

    img_end = ["\\end{center}",
               "\\end{figure}"]

    tbl_beg = ["\\begin{table}[htp]",
               "\\centering",
               "\\begin{tabular}{"+"*{"+str(columns+1)+"}{|>{\centering}p{"+str(15./(columns+1))[:3]+"cm}}|} ",
               "\\hline"]

    tbl_end = ["\\end{tabular}",
               "\\end{table}"]

    eqn_beg = ["\\begin{equation}"]
    
    eqn_end = ["\\end{equation}"," "]
    
    vid_beg = ["\\begin{figure}[ht]",
               "\includemovie[",
               "poster,"]
    
    vid_end = ["\end{figure}"]
    
    IMAGES,TABLES,VIDEOS,EQNS,DOTS="","","","",""

    qqq=0
    for i in media['image']:
        for j in img_beg:
            IMAGES=IMAGES+j+" \n "

        IMAGES = IMAGES+ "\epsfig{file="+i+",scale = 0.8} \n"

        for j in img_end:
            IMAGES=IMAGES+j+" \n "
        qqq=qqq+1
        if qqq%15==1:
            IMAGES = IMAGES + "\clearpage \n"

    for i in media['table']:
        for j in tbl_beg:
            TABLES=TABLES+j+" \n "

        tmp = ["& $ U_"+str(q+1)+ " $ " for q in range(self.mon)]
            
        TABLES = TABLES + " * "
        for k in tmp:
            TABLES = TABLES+k

        TABLES = TABLES + " \n \\tabularnewline \n \\hline \n "
        TABLES = TABLES + i +" \n "

        for j in tbl_end:
            TABLES=TABLES+j+" \n "

    for i in media['dot']:
        MySys("$DOTPATH./ladot "+i+".ladot")
        DOTS=DOTS+"\\input{"+i+".tex} \n"
   
def PublishImageList(name,title,filelist):
    f=open(filelist)
    media={'table': [], 'image': [], 'video': [], 'dot': [], 'eqn': []}
    for line in f:
        media['image'].append(line.replace('\n',''))
    Publish(name,title,media)

def MeshTraj(mesh,velocity,dt,tsteps,record=0,interactive=1):

    if record !=0:
        nm=cname(1)
        
    for j in range(tsteps):
        for x in mesh.coordinates(): 
            xc = TStep(velocity,x,dt)              
            for k in range(len(x)):
                x[k] = xc[k]
        if j%4500==0 and interactive==1:
            p=plot(mesh,interactive=True,axes=True)
        else:
            p=plot(mesh,axes=True)
        if record !=0:
            if j%25==0:
                p.write_png(record+"_"+nm.next()+".png")
            
def TStep(velocity,x,dt):
    try:
        a=velocity(x)
    except IndexError:
        a=[velocity[i](x) for i in range(len(velocity))]
    return [x[i]+dt*a[i] for i in range(len(x))]

def MakeInitialCondition(mesh,Q,MyFile,hist='u_2'):

    f = TFile(MyFile)
    u4 = gDirectory.Get(hist)
    u4.Draw()
    u0 = Function(Q)
    for i,j in enumerate(mesh.coordinates()[:]):
        tmp = j[0]*j[1]/j[2]
        tmp = u4.FindBin(tmp)
        u0.vector()[i]=u4.GetBinContent(tmp)

    return u0
        

