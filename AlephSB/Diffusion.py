
## @package Diffusion
# Class for solution of diffusion equations over Dynamical Systems in Dolfin.
#????????????????????????????????????????????????????????????????????????????
#
#
#  We seek to solve 
# \f$\partial_t u = \partial_i (D_{ij}\partial_j u - F_i u)\f$.
#  The functional form is found by multiplying a test function v and integrating over the mesh.
# \f$\int_\Omega v \partial_t u dx =\int_\Omega v  \partial_i (D_{ij}\partial_j u - F_i u) dx \f$
# Continuing with integration by parts yields \f$ \int_{\partial\Omega} v (D_{ij}\partial_j u - F_i u) n_i ds - \int_\Omega \partial_i v (D_{ij}\partial_j u - F_i u) dx\f$. 
# Here \f$D_{ij}\f$ denotes the spatial dependance of the Diffusion, \f$F_{i}\f$, is the drift force, \f$n_i\f$, a normal vector to the mesh.  The mesh is denoted \f$\Omega\f$ with boundary \f$\partial\Omega \f$.
#
# We approximate to find an expression for \f$u_{t+\delta t} \f$ given the solution at \f$u_t\f$.
#
# \f$\int_\Omega v (u_{t+\delta t} -  u_t) dx \approx \delta t [\int_{\partial\Omega} v (D_{ij}\partial_j u - F_i u) n_i ds - \int_\Omega \partial_i v (D_{ij}\partial_j u - F_i u) dx ]\f$
# We take the average of the expression [...] over time \f$\delta t\f$ by asserting \f$u \approx \frac{u_{t+\delta t} +  u_t}{2}\f$.  Continuing gives 
#
# \f$\frac{\delta t}{2} [\int_{\partial\Omega} v (D_{ij}\partial_j u_{t+\delta t} - F_i u_{t+\delta t}) n_i ds - \int_\Omega \partial_i v (D_{ij}\partial_j u_{t+\delta t} - F_i u_{t+\delta t}) dx + \f$ \f$ \int_{\partial\Omega} v (D_{ij}\partial_j u_t - F_i u_t) n_i ds - \int_\Omega \partial_i v (D_{ij}\partial_j u_t - F_i u_t) dx ]\f$
#
# Collecting terms with \f$u_{t+\delta t} \f$ on one side of equality and \f$u_t\f$ on the other gives the variational forms a and L that we will use.
#
# \f$ a = \int_\Omega (v u_{t+\delta t}+\frac{\delta t}{2}\partial_i v [D_{ij}\partial_j u_{t+\delta t} - F_i u_{t+\delta t}]) dx \f$ \f$- \frac{\delta t}{2} \int_{\partial\Omega} v (D_{ij}\partial_j u_{t+\delta t} - F_i u_{t+\delta t}) n_i ds   \f$ 
#
# \f$ L = \int_\Omega (v u_{t}-\frac{\delta t}{2}\partial_i v [D_{ij}\partial_j u_{t} - F_i u_{t}]) dx + \frac{\delta t}{2} \int_{\partial\Omega} v (D_{ij}\partial_j u_{t} - F_i u_{t}) n_i ds   \f$.
#
# Compare with the lines of python code:
#
# \f$a = v*u*dx + 0.5*k*(inner(grad(v),velocity*u)*dx + inner(grad(v), D*grad(u))*dx -v*inner(D*grad(u)-velocity*u,n)*ds)\f$
#
# \f$L = v*u0*dx - 0.5*k*(inner(grad(v),velocity*u0)*dx + inner(grad(v), D*grad(u0))*dx-v*inner(D*grad(u0)-velocity*u0,n)*ds)\f$
#
#@author Nathan Borggren
#
#????????????????????????????????????????????????????????????????????????????
from Dynamics import *
from Noise import *
from numpy import abs as ABS
from numpy import log as LN
from random import choice
from Dissonance import *
from Analysis import cname

## The diffusion class for study of \f$\partial_t u = \partial_i (D_{ij}\partial_j u - F_i u)\f$. 


class DirichletBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary

def boundary(x):
    return x[0] < DOLFIN_EPS-1 or x[1] < DOLFIN_EPS-1 

def bc_0(x):
    return pow(pow(x[0]-0.14175347706,2)+pow(x[1]-462.572738749,2),0.5) < 10

def bc_1(x):
    return pow(pow(x[0]-644.051823392,2)+pow(x[1]-0.00611657521734,2),0.5) < 10

def bc_2(x):
    return pow(pow(x[0]-29.9798485002,2)+pow(x[1]-178.229919659,2),0.5) < 10

#class DirichletBoundary(SubDomain):
#    def inside(self, x, on_boundary):
#        return on_boundary
class Diffusion(Dynamics, Noise, Dissonance):
    
    ## Initialize Defaults
    def Init(self):

        Dynamics.Init(self)
        Noise.Init(self)
        self.normT = array('d')
        self.norm = array('d')
        self.sMesh=[]

        ##sets the necessary variables self.mesh and self.Q from a file input or generates two d from set ranges.   
    def Mesh(self,mesh_file,rewrite=False,dim = 1,vdim=2):

        if rewrite==True:
            MySys("rm "+mesh_file)
        if self.Log==1:
            MySys("echo mesh used is named: "+mesh_file+" >> "+self.output+".log")
            MySys("echo mesh started at: " +strftime('%y.%m.%d -- %H.%M.%S ')+" >> "+self.output+".log")
        try:
            mesh = Mesh(mesh_file)
        except RuntimeError:
            mesh = UnitSquare(self.nx,self.nx)
                            
            for k in mesh.coordinates():
                k[0] = (self.maxs[0]-self.mins[0])*k[0]+self.mins[0]
                k[1] = (self.maxs[1]-self.mins[1])*k[1]+self.mins[1]

            mesh_out = File(mesh_file)
            mesh_out << mesh
            print "generating the mesh then"

        if self.Log==1:
            MySys("echo mesh loaded at: " +strftime('%y.%m.%d -- %H.%M.%S ')+" >> "+self.output+".log")


        #qRectangle(100,100,.97,.97,mins=[0,0],maxs=[1000,1000])
        #self.mesh = mesh
        #self.Q = FunctionSpace(mesh, "Lagrange", 2)
        #self.V = VectorFunctionSpace(mesh, "Lagrange", 2)
        return

    ##solves the steady state equation for given parameters and rights to output and screen. (no it does not)
    def Steady(self,output,t0='hist'):
        if self.Log==1:
            MySys("echo Simulation begun at: " +strftime('%y.%m.%d -- %H.%M.%S ')+" >> "+self.output+".log")

        mesh = self.mesh
        velocity = self.velocity

        D,velocity = self.dpar['D'], self.velocity
        g  = Constant(0.0)
        bc = DirichletBC(self.Q, g, boundary)

        # Create FunctionSpaces
        Q = FunctionSpace(mesh, "CG", 1)
        V = VectorFunctionSpace(mesh, "CG", 2)
        n = FacetNormal(mesh)

        # Initialise source function and previous solution function
        u0 = self.GetInit(Q,t0)
        u0 = self.Normalize(u0)
        du,v,u1,u,zero = TrialFunction(Q), TestFunction(Q),Function(Q),Function(Q),Function(Q)
        lam1,lam2 = Constant(5000000.),Constant(50000000.0)

        nn = len(u1.vector()[:])
        for i in range(nn):
            u1.vector()[i]=u0.vector()[i]

        D = 0.005 
        L = u*v*dx+inner(grad(v), D*grad(u)-velocity*u)*dx-v*inner(D*grad(u)-velocity*u,n)*ds- v*u1*dx
        a = derivative(L, u, du)
  
        A = assemble(a)

        b = assemble(L)
        bc.apply(A,b)        
        solve(A, u.vector(), b)

        #uu = self.Zoom(u,[0,3],[1,5.2])
        #        a = u*v*dx+inner(grad(v), D*grad(u)-velocity*u)*dx-v*inner(D*grad(u)-velocity*u,n)*ds
        #L = lam1*u*inner(-velocity,n)*ds+lam2*inner(D*grad(u)-velocity*u,n)*ds
        #L = lam1*u*dx
        #L=u*u*dx
        #L=lam2*inner(D*grad(v)-velocity*v,n)*ds+lam1*v*dx
        #       L= v*u0*dx 
        #      A = assemble(a)

        #     b = assemble(L)
        #    solve(A, u1.vector(), b)
        #pp = plot(uu)
        p=plot(u)
        interactive()
        
        p.write_png(output+"_Steady.png")
        out_file = File(output+"s.pvd")
        out_file << u1

        self.SaveFunction(u1,output+"_steady") 

        return

    ##runs a trajectory from a given initial probability distribution.  Defaults to initial condition from histogram object, e.g. filled with langevin trajectories.
    def Run(self,output,t0='expr',dim=1,vdim=2,viddy=True,ROOT=False):

        print self.nx
        self.mesh =  Rectangle(self.mins[0],self.mins[1],self.maxs[0],self.maxs[1],self.nx,self.nx) #why?
        meshtmp =   Rectangle(self.mins[0],self.mins[1],self.maxs[0],self.maxs[1],int(self.nx/3),int(self.nx/3)) #why?
        
        mesh=self.mesh
        ggg = File(output+"_mesh.xml")
        ggg<<mesh

        if self.MyLog==1:
            print self.output
            print strftime('%y.%m.%d -- %H.%M.%S ')
            MySys("echo Simulation begun at: " +strftime('%y.%m.%d -- %H.%M.%S ')+" >> "+self.output+".log")
            MySys("echo diffusion parameters, t0,tf,dt are "+str(self.dpar)+" >> "+self.output+".log")
            MySys("echo parameters are "+str(self.par)+" >> "+self.output+".log")

            MySys("echo mesh used is "+output+"_mesh.xml  >> "+self.output+".log")
        velocity = self.velocity
        D = self.D

        t,T,k = self.dpar['t'], self.dpar['T'],self.dpar['k']
        print self.dpar

        # Create FunctionSpaces
        Q = self.Q
        V = self.V

        # Set up boundary condition
        g  = Constant(0.00001)
        bc = DirichletBC(Q, g, boundary)
        #bc = DirichletBC(Q, g, DirichletBoundary.inside)

        # Initialise source function and previous solution function
        #print Q,self.nu0
        u0=self.u0
        #u0 = self.GetInit(Q,t0)
        #plot(u0,interactive=True)
        #u0 = self.Normalize(u0)
        u,v,u1,u2 = TrialFunction(Q), TestFunction(Q), Function(Q), Function(Q)

        # a = v*u*dx + 0.5*k*(inner(grad(v),-velocity*u)*dx + inner(grad(v), D*grad(u))*dx-v*inner(D*grad(u)-velocity*u,n)*ds)
        # L = v*u0*dx - 0.5*k*(inner(grad(v),-velocity*u0)*dx + inner(grad(v), D*grad(u0))*dx-v*inner(D*grad(u0)-velocity*u0,n)*ds)
        
        a = v*u*dx + 0.5*k*(inner(grad(v),-velocity*u)*dx + inner(grad(v), D*grad(u))*dx)
        L = v*u0*dx - 0.5*k*(inner(grad(v),-velocity*u0)*dx + inner(grad(v), D*grad(u0))*dx)

        A = assemble(a)

        j,jj=0,0
        counter=cname(1)

        #cN = TCanvas()
        #edensity=[]
        # Time-stepping
        while t < T:

            b = assemble(L)
            bc.apply(A,b)
            solve(A, u1.vector(), b)
            #u1.vector()[:] /= 
            print 'norm', assemble(u1*dx, mesh=self.mesh)
 
            #X = Function(Q,"x[0]") #for i in range(2)]
            #ops = [Expression("x[0]"),Expression("x[1]")]
            #for i,j in zip(fops,ops):
            #    i.interpolate(j)

            #observables = [dot(velocity,-D*grad(u1)+velocity*u1)/dot(-D*grad(u1)+velocity*u1,-D*grad(u1)+velocity*u1)*dx,X*dx,fops[1]*u1*dx]

            #print '<F.J>', assemble(dot(velocity,-D*grad(u1)+velocity*u1)*dot(-D*grad(u1)+velocity*u1,-D*grad(u1)+velocity*u1)*dx, mesh=self.mesh)
            #print assemble(o,mesh=self.mesh),"next" 

            #u1=self.Normalize(u1,t=t)
            # Copy solution from previous interval
            u0.assign(u1)

            if viddy == True:
                if ROOT==False:
                    p=plot(u1)
                    if jj%self.png==0:
                        j+=1
                    #use cname, fix this crap
                        tmp = (6-len(str(j)))*"0"+str(j)
                else:
                    pass
                    #edensity.append(energy_density(mesh,u0))
                    #tmpE=energy_density(meshtmp,u0)
                    #tmpE.Draw("surf2")
                #u2=self.LR(u1,u2)
                if self.smesh==1:
                   
                    uu.interpolate(u1)
                    
                    r=plot(uu)
                    r.write_png(self.output+"_Z.png")              


                #p=plot(u1)
                #q=plot(u2)

                #p.write_png(output+"_"+tmp+"_P.png")
                #q.write_png(output+"_"+tmp+"_S.png")              

            if jj%self.xml==0:
                tmp=counter.next()
                #self.SaveFunction(u1,output+"_"+tmp) 
                ufile = File("./"+output+"_"+tmp+".xml")
                ufile << u1
            t += k
            jj += 1

            #tmpE=Func2Graph(meshtmp,u0)
            #tmpE.Draw("surf2")

        zwom = input("sheat")
        #nn = len(self.normT)
        #self.renorm = TGraph(nn,self.normT,self.norm)

        if self.MyLog==1:
            MySys("echo everything short of "+output+"_"+counter.next()+".xml written  >> "+self.output+".log")
            MySys("tar -cvf "+output+".tar "+output+".log "+output+"_*")
            MySys("rm "+output+"_*.gz "+output+"_*.png ")
            MySys("echo Simulation finished at: " +strftime('%y.%m.%d -- %H.%M.%S ')+" >> "+self.output+".log")


    def FirstPassage(self,bc="bc_0"):

        #mesh=self.mesh
        mesh =  Rectangle(self.mins[0],self.mins[1],self.maxs[0],self.maxs[1],500,500)
        def bnd_123_0(x):
            #return pow(pow(x[0]-0.14175347706,2)+pow(x[1]-462.572738749,2),0.5) < 10
            return pow(pow(x[0]-0.14175347706,2)+pow(x[1]-362.572738749,2),0.5) < 10
        #meshtmp =   Rectangle(self.mins[0],self.mins[1],self.maxs[0],self.maxs[1],50,50) #why?
        u0 = Constant(0.0)
        velocity = self.velocity
        D = self.D
        #Q=self.Q
        Q=FunctionSpace(mesh,"Lagrange",2)
# Define variational problem


        v = TestFunction(Q)
        u = TrialFunction(Q)
        a = v*inner(velocity,grad(u))*dx-inner(grad(v), D*grad(u))*dx
        L = -v*dx
        bc_tmp = DirichletBC(Q, u0, bnd_123_0)
        problem = VariationalProblem(a, L, bc_tmp)
        u = problem.solve() 
        #plot(u,interactive=True)
        #tmpFPT=Func2Graph(mesh,u)
        #tmpFPT.Draw("surf2")

        for i in range(100):
            print u(10*i,.001)
        if self.MyLog==1:
            print self.output
            print strftime('%y.%m.%d -- %H.%M.%S ')
            MySys("echo Simulation begun at: " +strftime('%y.%m.%d -- %H.%M.%S ')+" >> "+self.output+".log")
            MySys("echo fpt from x to: "+str(X)+" "+str(Y)+" >> "+self.output+".log")
            MySys("echo parameters are "+str(self.par)+" >> "+self.output+".log")

            MySys("echo mesh used is "+self.output+"_mesh.xml  >> "+self.output+".log")

    ##zooms in on a particular cranny or nook and appends submesh to sMesh.
    def Zoom(self,xrng,yrng,nx=250):
    # Structure sub domain
        smesh = UnitSquare(nx,nx)

        dx = xrng[1]-xrng[0]
        dy = yrng[1]-yrng[0]

        for i in smesh.coordinates():
            i[0]=i[0]*dx+xrng[0]
            i[1]=i[1]*dy+yrng[0]

        self.sMesh.append(smesh)

        return 

    def SaveFunction(self,anyF,name):
        zz=File(name+".xml")
        zz<<anyF.vector()[:]
        MySys("gzip "+name+".xml")

    def Normalize(self,u,t=-1):
        norm = assemble(u*dx, mesh=self.mesh)
        u.vector()[:] /= norm
        if t>0:
            self.normT.append(t)
            self.norm.append(norm)
        return u

    def GetInit(self,Q,smooth = 5, expr = "exp(-pow(x[0]-200.5,2)-pow(x[1]-200.5,2))", func = None,n=5):
        print "getting init"
        u0 = Function(Q)
        u1 = Function(Q)

        if self.t0 == 'hist':
            u0 = self.hist2u0(u0)

        if self.t0 == 'traj':
            u0 = self.traj2u0(n,u0,u1)

        if self.t0 == 'func':
            try:
                z = Function(Q,self.nu0)
                u0 = self.func2u0(z,u0)
            except TypeError:
                self.u0 = self.nu0
                return self.u0
        if self.t0 == 'rand':
            u0,u1 = Function(Q),Function(Q)
            u0 = self.rand2u0(u0,u1,n=self.rand)

        if self.t0 == 'expr':
            utmp = Expression(expr)
            u0 = Function(Q)
            u0.interpolate(utmp)            

        return u0

    def func2u0(self,MyFunc,uu):
        uu.interpolate(MyFunc)
        return uu
    def expr2u0(self,Q,expr="exp(-pow(x[0]-2.5,2)-pow(x[1]-2.5,2))"):
        self.u0 = project(Expression(expr),Q)

    def hist2u0(self,u0):

        nn = len(self.mesh.coordinates()[:])
        for i in range(nn):
            tmp=self.mesh.coordinates()[i]
            mybin = self.hist.FindBin(tmp[0],tmp[1])
            val = self.hist.GetBinContent(mybin)
            u0.vector()[i]=val
        return u0

    def traj2u0(self,n,u0,u1):
        tt,xx = self.Traj(self.rzero(),fast=0)
        nn=len(tt)
        for i in range(n):
            tmp = choice(range(nn))
            a1, b1 = max([xx[0][tmp],0.001]),max([xx[1][tmp],0.001])
            it = dict(a = a1,b=b1)
            gaus = Expression("exp(-(pow(x[0]-a,2)+pow(x[1]-b,2))/0.001)",defaults=it)
            u1.interpolate(gaus)
            for j in range(len(u0.vector()[:])):
                u0.vector()[j]+=u1.vector()[j]
        return u0

    def rand2u0(self,u0,u1,n=7):
        nn=len(u0.vector())
        
        for i in range(n):
            xx=self.rzero()
            it = dict(a = xx[0],b=xx[1])
            #gaus = Expression("exp(-(pow(x[0]-a,2)+pow(x[1]-b,2))/0.01)",defaults=it)
            gaus = Expression("exp(-(pow(x[0]-a,2)+pow(x[1]-b,2))/0.001)",defaults=it)
            u1.interpolate(gaus)
            for j in range(len(u0.vector()[:])):
                u0.vector()[j]+=u1.vector()[j]

        return u0

    #def GetCalc(self):
    #    self.MyCalc()
    #    self.mainloop()
              
    def Clean(self,output,dest):
        MySys("tar -czvf "+ output + "_dat.tar.gz " + output + ".pvd " + output+"*.vtu")
        MySys("mv "+ output + "_dat.tar.gz "+ output + ".pvd " + output+"*.vtu " + dest)
        #MySys("rm "+ self.name+".pvd " + self.name+"*.vtu")

        MySys("tar -czvf "+ output + "_P_png.tar.gz " + output + "*_P.png")
        MySys("rm "+ output + "_P.mp4")
        MySys("ffmpeg -qscale 5 -r 5 -b 9600 -i "+output+"_%06d_P.png "+output+"_P.mp4")
        MySys("mv "+ output + "_P_png.tar.gz "+ output + "*_P.png " + dest)
        #MySys("rm "+ output + "*.png")

        MySys("tar -czvf "+ output + "_S_png.tar.gz " + output + "*_S.png")
        MySys("rm "+ output + "_S.mp4")
        MySys("ffmpeg -qscale 5 -r 2 -b 9600 -i "+output+"_%06d_S.png "+output+"_S.mp4")
        MySys("mv "+ output + "_S_png.tar.gz "+ output + "*_S.png " + dest)

        MySys("tar -czvf "+ output + "_Z_png.tar.gz " + output + "*_Z.png")
        MySys("rm "+ output + "_Z.mp4")
        MySys("ffmpeg -qscale 5 -r 5 -b 9600 -i "+output+"_%06d_Z.png "+output+"_S.mp4")
        MySys("mv "+ output + "_Z_png.tar.gz "+ output + "*_Z.png " + dest)

        MySys("tar -czvf "+ output + "_xml.tar.gz " + output + "*.xml")
        MySys("mv "+ output + "_xml.tar.gz "+ output + "*.xml " + dest)
        #MySys("rm "+ output + "*.xml")

        if self.Log==1:
            MySys("echo wrote "+output+ " files to dir: "+dest+" >> "+self.output+".log")
        return
        

