##@package Dynamics
#Classes for Dynamical Systems 
##@author Nathan Borggren

from dolfin import *
from ROOT import TH3F,TH2F,TH1F, TF1, TGraph, TGraph2D, gStyle, TPolyLine3D, TCanvas, TF3, TPolyLine, TPostScript
from array import array as MyArray
from time import time, strftime, gmtime
from Algebra import *
from random import random
from numpy.fft import *
from Analysis import *
from Calc import *

gStyle.SetPalette(1)  

##The Dynamical system class
class Dynamics(Algebra,Calc):
    def Init(self,S=False,Ti=[0.,15.,50000],nx=256,
             mins = [-30,-30,0],
             maxs = [30,30,50],
             doms = [[0.2,1.8],[-10,4.5],[-2.5,6],[-50,150],[-3.5,15]],
             ROOT = 1,dof=3,mon=5,system='Lorenz',
             f=["-sigma*x[0]+sigma*x[1]",
                "rho*x[0]-x[1]-x[0]*x[2]",
                "-beta*x[2]+x[0]*x[1]"],
             u=["1","y/x","x/y","x*z/y","x*y/z"],
             par={'sigma':10,'rho':28.,'beta':8./3.},
             spar = [symbol("\\sigma"),symbol("\\rho"),symbol("\\beta")],
             A=matrix([[10.,10.,0,0,0],
                       [-1,0,28.,-1,0],
                       [-8./3.,0,0,0,1]]),
             B=matrix([[0,0,0],
                       [-1,1,0],
                       [1,-1,0],
                       [1,-1,1],
                       [1,1,-1]]),
             imagine=0,
             hists=[],
             MyLog=1):

        self.u = u
        self.nx=nx
        self.mins=mins
        self.maxs=maxs
        self.doms=doms
        self.ROOT=ROOT
        self.imagine=imagine

        Algebra.Init(self,S=S,dof=dof,mon=mon,system=system,f=f,u=u,par=par,spar=spar,A=A,B=B,imagine=imagine)        

        # if ROOT == 0:
        #     try:
        #         self.velocity = Expression((self.f[i] for i in range(self.dof)),defaults=self.par)
        #     except TypeError:
        #         self.GetVelocity()
            
        self.npSect=0
        self.pSect=[]
        #self.T(Ti[0],Ti[1],Ti[2])
        self.nLines = []

    ##Flow field calculated at x.
    def Flow(self,x):
        if self.ROOT == 0:
            return [self.velocity(x)[i] for i in range(self.dof)]
        elif self.S == False:
            sigma=self.par['sigma']
            rho=self.par['rho']
            beta=self.par['beta']
            return [-sigma*x[0]+sigma*x[1],rho*x[0]-x[1]-x[0]*x[2],-beta*x[2]+x[0]*x[1]]
        else:
            sigma=self.spar[0]
            rho=self.spar[1]
            beta=self.spar[2]
            
            return [-sigma*x[0]+sigma*x[1],rho*x[0]-x[1]-x[0]*x[2],-beta*x[2]+x[0]*x[1]]              

    ##Set time specifications for trajectories
    def T(self,ti,tf,tsteps):
        self.ti=ti
        self.tf=tf
        self.dt=(tf-ti)/float(tsteps)
        self.tsteps=tsteps
        self.mons(l=[ti,tf,tsteps])
        self.space(l=[ti,tf,tsteps])
 
    def mons(self,l=[0,15,50000]):
        if self.imagine==1:
            self.hUr = [TH1F("u_"+str(k)+"_r","u_"+str(k)+"_r",2500,self.doms[k][0],self.doms[k][1]) for k in range(self.mon)]
            self.tUr = [TH1F("tu_"+str(k)+"_r","tu_"+str(k)+"_r",int(l[2]),l[0],l[1]) for k in range(self.mon)]

            self.hUi = [TH1F("u_"+str(k)+"_i","u_"+str(k)+"_i",2500,self.doms[k][0],self.doms[k][1]) for k in range(self.mon)]
            self.tUi = [TH1F("tu_"+str(k)+"_i","tu_"+str(k)+"_i",int(l[2]),l[0],l[1]) for k in range(self.mon)]

        else:
            self.hU = [TH1F("u_"+str(k),"u_{"+str(k)+"}",2500,self.doms[k][0],self.doms[k][1]) for k in range(self.mon)]
            #print int(,l[0],l[1]
            self.tU = [TH1F("tu_"+str(k),"tu_{"+str(k)+"}",int(l[2]),l[0],l[1]) for k in range(self.mon)]

            if self.dof==3:
                self.hist = TH3F("histr","histr",100,self.mins[0],self.maxs[0],100,self.mins[1],self.maxs[1],100,self.mins[2],self.maxs[2])
 
    def space(self,l=[0,15,50000]):      
        
        if self.imagine==1:
            self.hxr = [TH1F("x_"+str(k)+"_r","Re(x_{"+str(k)+"})",2500,self.mins[k],self.maxs[k]) for k in range(self.dof)]
            self.txr = [TH1F("tx_"+str(k)+"_r","Re(tx_{"+str(k)+"})",int(l[2]),l[0],l[1]) for k in range(self.dof)]
            self.hxi = [TH1F("x_"+str(k)+"_i","Im(x_{"+str(k)+"})",2500,self.mins[k],self.maxs[k]) for k in range(self.dof)]
            self.txi = [TH1F("tx_"+str(k)+"_i","Im(tx_{"+str(k)+"})",int(l[2]),l[0],l[1]) for k in range(self.dof)]
            if self.dof==3:
                self.histr = TH3F("histr","histr",100,self.mins[0],self.maxs[0],100,self.mins[1],self.maxs[1],100,self.mins[2],self.maxs[2])
                self.histi = TH3F("histi","histi",100,-self.maxs[0],self.maxs[0],100,-self.maxs[1],self.maxs[1],100,-self.maxs[2],self.maxs[2])

            if self.dof==2:
                self.histr = TH2F("histr","histr",100,self.mins[0],self.maxs[0],100,self.mins[1],self.maxs[1])
                self.histi = TH2F("histi","histi",100,-self.maxs[0],self.maxs[0],100,-self.maxs[1],self.maxs[1])


        else:
            self.hx = [TH1F("x_"+str(k),"x_{"+str(k)+"}",100,self.mins[k],self.maxs[k]) for k in range(self.dof)]
            self.tx = [TH1F("tx_"+str(k),"tx_{"+str(k)+"}",int(l[2]),l[0],l[1]) for k in range(self.dof)]

            if self.dof==3:
                
                self.hist = TH3F("hist","hist",2500,self.mins[0],self.maxs[0],2500,self.mins[1],self.maxs[1],2500,self.mins[2],self.maxs[2])
 

    ##Time iteration, defaults to runge kutte
    def Tstep(self, x, dt, fast=0, q = 1, t = 0):
        #print time stepping
        if fast == 1:
            a = self.Flow(x)
            #print a, 'a'
            return [x[i]+(q*t-t+dt)*a[i] for i in range(self.dof)]

        elif q==0 and fast == 0:
            x1 = x
            a1 = self.Flow(x1)
        
            x2 = [x[i] + 0.5*a1[i]*dt for i in range(self.dof)]
            a2 = self.Flow(x2)
        
            x3 = [x[i] + 0.5*a2[i]*dt for i in range(self.dof)]
            a3 = self.Flow(x3)
        
            x4 = [x[i] + 0.5*a3[i]*dt for i in range(self.dof)]
            a4 = self.Flow(x4)
 
            return [x[i] + (dt/6.0)*(a1[i]+2*a2[i]+2*a3[i]+a4[i]) for i in range(self.dof)]

    ## Returns arrays of trajectories and perhpas monomials for initial condition x 
    def Traj(self,x,Fill=None,U=1,Plot=0,fast=1,store=0,q=0,gen=False,noise=0):

        print "Im trying 1" 
        if self.imagine==1:
            trajr=[array('d') for i in range(self.dof)]
            traji=[array('d') for i in range(self.dof)]

        else:
            traj=[array('d') for i in range(self.dof)]

        xraj=array('d')   
        xraj.append(self.ti)

        if self.imagine==0:
            ost = ""
            for j in range(self.dof):
                ost = ost+"x["+str(j)+"],"
            ost = ost[:-1]
        else:
            ostr, osti = "",""
            for j in range(self.dof):
                ostr,osti = ostr+"x["+str(j)+"].real,",osti+"x["+str(j)+"].imag,"
                
            ostr, osti = ostr[:-1],osti[:-1]

        print "Im trying 2"
        if U==1:

            tmp=self.Mons(x)

            if self.imagine==1:
                utrajr=[array('d') for i in range(self.mon)]
                utraji=[array('d') for i in range(self.mon)]

                for i in range(self.mon):
                    utrajr[i].append(tmp[i].real)
                    utraji[i].append(tmp[i].imag)

            else:
                utraj=[array('d') for i in range(self.mon)]
                for i in range(self.mon):
                    utraj[i].append(tmp[i])

        for i in range(self.dof):
            if self.imagine==1:
                trajr[i].append(x[i].real)
                traji[i].append(x[i].imag)
            else:    
                traj[i].append(x[i])

        if gen==True:
            z=qIntegrate(x,self.dt,self.velocity)
       
        for j in range(self.tsteps):
            if gen == False:
                if noise==0:
                    x = self.Tstep(x,self.dt,fast,q=q)
                else:
                    xtmp = self.Tstep(x,self.dt,fast,q=q)
                    x = self.NoiseStep(xtmp,self.dt)
            else:
                x = z.next()
                
            xraj.append((j+1)*self.dt+self.ti)
            for i in range(self.dof):
                if self.imagine==1:
                    trajr[i].append(x[i].real)
                    traji[i].append(x[i].imag)
                else:    
                    traj[i].append(x[i])

            if U==1:
                tmp=self.Mons(x)
                for i in range(self.mon):
                    if self.imagine==1:
                        utrajr[i].append(tmp[i].real)
                        utraji[i].append(tmp[i].imag)
                        self.hUr[i].Fill(tmp[i].real)
                        self.hUi[i].Fill(tmp[i].imag)
                        #self.tUr[i].Fill((j+1)*self.dt+self.ti,tmp[i].real)
                        self.tUr[i].AddBinContent(j+1,tmp[i].real)
                        #self.tUi[i].Fill((j+1)*self.dt+self.ti,tmp[i].imag)
                        self.tUi[i].AddBinContent(j+1,tmp[i].imag)
                    else:
                        utraj[i].append(tmp[i])
                        self.hU[i].Fill(tmp[i])
                        self.tU[i].Fill((j+1)*self.dt+self.ti,tmp[i])

            for i in range(self.dof):

                if self.imagine==0:
                    self.hx[i].Fill(x[i])
                    #self.tx[i].Fill((j+1)*self.dt+self.ti,x[i])
                    self.tx[i].AddBinContent(j+1,x[i])
                    exec("self.hist.Fill("+ost+")")
                else:
                    self.hxr[i].Fill(x[i].real)
                    #self.txr[i].Fill((j+1)*self.dt+self.ti,x[i].real)
                    self.txr[i].AddBinContent(j+1,x[i].real)
                    #self.txr[i].Fill(x[i].real)
                    self.hxi[i].Fill(x[i].imag)
                    #self.txi[i].Fill((j+1)*self.dt+self.ti,x[i].imag)
                    self.txi[i].AddBinContent(j+1,x[i].imag)
                    #self.txi[i].Fill(x[i].imag)
                    exec("self.histr.Fill("+ostr+")")
                    exec("self.histi.Fill("+osti+")")


            # for i in range(self.npSect):
            #     print self.npSect, 'npSect'
            #     distance = self.dist(self.surf[i],x)
            #     if distance < .1:
            #         fill = self.fill(i,x)
            #         exec("self.pSect[i].Fill("+str(fill[0])+','+str(fill[1])+")")
        
        self.last=x

        if Plot==1:
            if U==1:
                if self.imagine==1:
                    c1,c2,c4,c5=self.Prettify(xraj,[trajr,traji],uu=[utrajr,utraji])
                    c1.Draw()
                    c2.Draw()
                    c5.Draw()
                    c4.Draw()
                else:
                    c1,c2,c4=self.Prettify(xraj,traj,uu=utraj)
            else:
                c1,c2,c4=self.Prettify(xraj,traj)
            c1.Draw()
            #c2.Draw()
            #c5.Draw()
            c4.Draw()
            zwom=input("sheat")

        if U==1:
            try:
                return xraj, traj, utraj
            except UnboundLocalError:
                return xraj,[trajr,traji],[utrajr,utraji]
        else:
            return xraj, traj

    ## run multiple trajectories
    def nTraj(self,n,fast=1,U=0,gen=False,imagine=0,Plot=1,RunList=[],noise=0):
        try:
            self.hist
        except AttributeError:
            self.hReset()

        t=0
        before=time()

        if self.MyLog==1:
            MySys("echo "+strftime("%H_%M_%S")+" :Begin Trajectory >> "+self.output+".log")

        #x=self.rzero()
        #xx,yy,uu = self.Traj(x,Fill=1,U=1,fast=fast)
        #self.hU = [TH1F("U_"+str(i)) for i in range(self.mon)]
      
        for i in range(n):
            if i%10==0:
                t+=(time()-before)/60.
                before = time()
                print t, " minutes for ", i, "trajectories"

            if imagine==0:
                x=self.rzero()
                self.Traj(x,Fill=1,fast=fast,U=U,gen=gen,noise=noise)
            if imagine==1:
                x=self.rzero()
                tt,xx,uu = self.Traj(x,Fill=1,fast=fast,U=U,noise=noise)
            #if imagine==2:
            #    x=self.izero(1./2048.)
            #    tt,xx,uu = self.Traj(x,Fill=1,fast=fast,U=U)

                # if i%10==0:
                #     self.nLines.append(TPolyLine3D(len(tt),xx[0][0],xx[0][1],xx[0][2]))
                #     self.nLines[-1].SetLineColor(4)
                #     self.nLines.append(TPolyLine3D(len(tt),xx[1][0],xx[1][1],xx[1][2]))
                #     self.nLines[-1].SetLineColor(2)

        for x in RunList:
            print x
            tt,xx,uu = self.Traj(x,Fill=1,fast=fast,U=U)

        if self.MyLog==1:
            MySys("echo "+strftime("%H_%M_%S")+" :End Trajectory >> "+self.output+".log")            

        if U==1:
            if imagine==1:
                c1 = TCanvas()
                c1.Divide(3,2)
                c3 = TCanvas()
                c3.Divide(3,2)

                c2 = TCanvas()
                c2.Divide(2,2)
                c4 = TCanvas()
                c4.Divide(2,2)


                for i in range(self.mon):
                    c1.cd(i+1)
                    self.hUr[i].SetLineColor(4)
                    self.hUr[i].Draw()
                    self.hUi[i].SetLineColor(2)
                    self.hUi[i].Draw("same")
                    c3.cd(i+1)
                    self.tUr[i].SetLineColor(4)
                    self.tUr[i].Draw()
                    self.tUi[i].SetLineColor(2)
                    self.tUi[i].Draw("same")
                    
                for i in range(self.dof):
                    c2.cd(i+1)
                    self.hxr[i].SetLineColor(4)
                    self.hxr[i].Draw()
                    self.hxi[i].SetLineColor(2)
                    self.hxi[i].Draw("same")
                    c4.cd(i+1)
                    self.txr[i].SetLineColor(4)
                    self.txr[i].Draw()
                    self.txi[i].SetLineColor(2)
                    self.txi[i].Draw("same")           

            for a,b in enumerate(self.nLines):
                if a%2==0:
                    b.SetLineColor(4)
                else:
                    b.SetLineColor(2)
                # if a==0:
                #     b.Draw("APES")
                # else:
                #     b.Draw("PS")

            if self.imagine==1:
                c5 = TCanvas()
                c5.Divide(2)
                c5.cd(1)
                self.histr.Draw()
                c5.cd(2)
                self.histi.Draw()

                
                #c5.Print(self.name+"zoom_"+self.current.next()+".ps")

            if Plot==1:
                zwom=input("sheat")
            
        return 

    def FM_SYNTH(self,graph=None):
        pass

    def MeshTraj(self,tsteps,q,h,t=1,boundary=1,sz=12,smooth=0,inspect=1,shift=0):
      
        for j in range(tsteps):
            #print j,t
    
            mesh = self.mesh
            for x in mesh.coordinates():
                #print x
                xc = self.Tstep(x,h,fast=1) 
                #print xc
                for k in range(self.dof):
                    x[k] = xc[k]
                u = self.Mons(x)
                for i in range(self.mon):
                    if i<self.dof:
                        self.hx[i].Fill(x[i])
                    self.hU[i].Fill(u[i])
            t=q*t+h
            if smooth>0:
                self.mesh.smooth(smooth)
            if inspect==1 and j%200==0:
                pass
            if j%1000==1001:
                p=plot(mesh,interactive=True)
            else:
                p=plot(mesh,axes=True)

    ## Poincare sections to inspect
    def pSect(self,hSect,cSect,fSect):
        self.pSect=len(hSect)                 
        self.cSect=cSect  ## conditions to be satisfied for inclusion in histogram
        self.hSect=[TH2F("h"+str(i),"h"+str(i),self.nx,hSect[i][0],hSect[i][1]) for i in range(pSect)]
        self.fSect = fSect

    ## Reset histogram
    def hReset(self):
        if self.dof==2:
            self.hist = TH2F(self.name,self.name,self.nx,self.mins[0],self.maxs[0],self.nx,self.mins[1],self.maxs[1])
        if self.dof==3:
            self.histr = TH3F(self.name+'r',self.name+'r',self.nx,self.mins[0],self.maxs[0],self.nx,self.mins[1],self.maxs[1],self.nx,self.mins[2],self.maxs[2])
            self.histi = TH3F(self.name+'i',self.name+'i',self.nx,-self.maxs[0],self.maxs[0],self.nx,-self.maxs[1],self.maxs[1],self.nx,-self.maxs[2],self.maxs[2])

    ## reset flow field if parameters changed 
    def SetVelocity(self):
        if self.dof==2:
            self.velocity = Expression((self.f[0],self.f[1]),defaults=self.par)
        if self.dof==3:
            self.velocity = Expression((self.f[0],self.f[1],self.f[2]),defaults=self.par)

    ##how about a report
    def reset(self):
        self.SetDiffusion()
        #self.SetDimerization()
        #self.SetPartition()
        self.SetVelocity()
        self.u0 = Function(self.Q)
        #print self.su0.get()
        self.u0 = project(Expression(self.su0.get()),self.Q)
        
    def rzero(self,force=''):
        if force=='real':
            tmp = [(self.maxs[i]-self.mins[i])*random()+self.mins[i] for i in range(self.dof)] 
            return tmp
        if self.imagine==0:
            tmp = [(self.maxs[i]-self.mins[i])*random()+self.mins[i] for i in range(self.dof)]
        elif self.imagine==1:
            tmp = [(1+random()*1j)*(self.maxs[i]-self.mins[i])*random()+self.mins[i] for i in range(self.dof)]

        #tmp = [(1+random()*1j)*(self.maxs[i]-self.mins[i])*random()+self.mins[i] for i in range(self.dof)] #override
        self.logbook['inits'].append(tmp)
        return tmp

    def izero(self,width):
        tmp = [(self.maxs[i]-self.mins[i])*random()+self.mins[i] for i in range(self.dof)]
        tmpi = [(random()-1/2.)*width*1j +tmp[i] for i in range(3)]
        return tmpi   

    def SetSpaces(self):
        self.mesh = Rectangle(self.mins[0],self.mins[1],self.maxs[0],self.maxs[1],self.nx,self.nx)
        self.Q = FunctionSpace(self.mesh, "Lagrange", 2)
        self.V = VectorFunctionSpace(self.mesh, "Lagrange", 2)        
        self.space()


    def dzero(self):
        tmp_d0 = self.rzero()
        dtmp = -1
        #self.logbook.append(str(tmp))
        while self.refresh<1:
            dtmp = dtmp+1
            tmp_d0 = [tmp_d0[0],tmp_d0[1]*sin(2*pi*dtmp/8.),tmp_d0[2]*cos(2*pi*dtmp/8.)]
            yield tmp_d0
   
    def myFFT(self,x):
        tt,xx,uu = self.Traj(x,U=1,fast=0)
        n,gn=len(tt),len(uu)
        grsU = [TGraph(n,tt,uu[i]) for i in range(5)]         
        ffts = [fft(uu[i]).real for i in range(5)]
        ffts2= [fft(uu[i]).imag for i in range(5)]
        fftUt = [array('d') for i in range(5)]
        fftUt2 = [array('d') for i in range(5)]
        for i in range(5):
            for j in range(n):
                fftUt[i].append(ffts[i][j]**2+ffts2[i][j]**2)
                fftUt2[i].append(ffts[i][j])
        fftUr = [TGraph(n,tt,fftUt[i]) for i in range(5)]         
        fftUi = [TGraph(n,tt,fftUt2[i]) for i in range(5)]         
        c1 = TCanvas()
        c1.Divide(3,2)
        for i in range(5):
            c1.cd(i+1)
            fftUi[i].Draw("APE")
 
        fftsx1 = [fft(xx[i]).real for i in range(3)]
        fftsx2= [fft(xx[i]).imag for i in range(3)]
        fftsxA = [array('d') for i in range(3)]
        for i in range(3):
            for j in range(n):
                fftsxA[i].append(fftsx1[i][j])
        fftxr = [TGraph(n,tt,fftsxA[i]) for i in range(3)]         
        
        c2 = TCanvas()
        c2.Divide(3)
        
        for i in range(3):
            c2.cd(i+1)
            fftxr[i].Draw("APE")
        zwom = input("numbers continue")

    ## Make plots
    def Prettify(self,tt,xx,uu=None,store=0,more=0,hold=1,verbose=0):

        n=len(tt)


        if self.imagine == 0:
            grsX = [TGraph(n,tt,xx[i]) for i in range(self.dof)]

        else:
            grsX = [TGraph(n,tt,xx[0][i]) for i in range(self.dof)]
            grsXi = [TGraph(n,tt,xx[1][i]) for i in range(self.dof)]

        self.current.next()
        c2 = TCanvas()
        if self.dof==2:
            c2.Divide(2)
        else:
            c2.Divide(2,2)

        gColors=[2,4,6,7,8,9,11,12,13,14,15,16,17,18]
        gShapes=[27,26,23,22,28,15,11]

        for i in range(self.dof):
            c2.cd(i+1)
            grsX[i].SetLineColor(4)
            grsX[i].SetLineWidth(1)
            grsX[i].SetMarkerColor(4)
            grsX[i].SetMarkerStyle(7)
            grsX[i].SetTitle("x_{"+str(i)+"}")
            grsX[i].GetYaxis().SetTitle("x_{"+str(i)+"}")
            grsX[i].GetXaxis().SetTitle("time")
            grsX[i].Draw("APES")
            grsX[i].SetName("gX"+str(i)+self.current.next())
            self.root.append(grsX[i])

            if self.imagine==1:
                grsXi[i].SetLineColor(2)
                grsXi[i].SetLineWidth(1)
                grsXi[i].SetMarkerColor(2)
                grsXi[i].SetMarkerStyle(7)
                #grsXi[i].SetTitle("x_{"+str(i)+"}")
                #grsXi[i].GetYaxis().SetTitle("x_{"+str(i)+"}")
                #grsXi[i].GetXaxis().SetTitle("time")
                grsXi[i].Draw("PS")
                grsX[i].SetName("gX"+str(i)+self.current.next())
                self.root.append(grsXi[i])

        if verbose==1:c2.Print(self.dtime+self.name+"x_"+self.current.next()+".ps")

        if self.dof==3:
            if self.imagine==0:
                yy = TPolyLine3D(n,xx[0],xx[1],xx[2])
            else:
                yy = TPolyLine3D(n,xx[0][0],xx[0][1],xx[0][2])
                yyi = TPolyLine3D(n,xx[1][0],xx[1][1],xx[1][2])

        elif self.dof==2:
            yy = TGraph(n,xx[0],xx[1])
            print 'yy is ', yy,len(xx[0]),n
        elif self.dof==1:
            yy = TGraph(n,tt,xx[0])


        c4 = TCanvas()

        yy.SetLineColor(4)
        yy.Draw("APES")
        #yy.SetName("t"+self.current.next())

        if self.imagine==1:
            yyi.SetLineColor(2)
            yyi.Draw("PS")
            self.root.append(yyi)
            #yy.SetName("ti"+self.current.next())

            #zwom = input("sheat")
        if verbose==1:c4.Print(self.dtime+self.name+"trace_"+self.current.next()+".ps")
        self.media['image'].append(self.dtime+self.name+"trace_"+self.current.next()+".ps")
        self.root.append(yy)


        if uu!=None:

            if self.imagine==0:
                gn = len(uu)
                grsU = [TGraph(n,tt,uu[i]) for i in range(gn)]
            else:
                gn = len(uu[0])
                grsU = [TGraph(n,tt,uu[0][i]) for i in range(gn)]
                grsUi = [TGraph(n,tt,uu[1][i]) for i in range(gn)]

            c1 = TCanvas()
            c1.Divide(self.mon/2+1,2)

            for j,i in enumerate(grsU):
                i.SetLineColor(4)
                i.SetLineWidth(2)
                i.SetMarkerColor(4)
                i.SetName("gU"+str(j)+self.current.next())


            if self.imagine==1:
                for i in grsUi:
                    i.SetLineColor(2)
                    i.SetLineWidth(2)
                    i.SetMarkerColor(2)
                    i.SetName("gUi"+str(j)+self.current.next()) 

            for i in range(gn):
                c1.cd(i+1)
                grsU[i].SetTitle(self.u[i])
                grsU[i].GetYaxis().SetTitle(self.u[i])
                grsU[i].GetXaxis().SetTitle("time")
                self.root.append(grsU[i])
                grsU[i].Draw("APES")

                if self.imagine==1:
                    grsUi[i].Draw("PS")
                    self.root.append(grsUi[i])

            if verbose==1:c1.Print(self.name+"u_"+self.current.next()+".ps")

            if self.imagine==1:
                c5 = TCanvas()
                c5.Divide(2)
                c5.cd(1)
                self.histr.Draw()
                c5.cd(2)
                self.histi.Draw()
                if verbose==1:c5.Print(self.name+"zoom_"+self.current.next()+".ps")
                return c1, c2, c4, c5

            else:
                return c1,c2,c4

        return c2,c4
        if verbose==1:
            MySys("echo "+strftime("%H_%M_%S")+" :End Trajectory >> "+self.output+".log")            

    def sZ(self):

        Z=""
        if len(self.u)<=len(self.par):
            for i,j in enumerate(self.u):
                Z=Z+"g"+str(i)+"*"+j+"+"
            Z=Z[:-1]
            print Z
            self.partition=Z
            return Expression(Z,defaults=self.par)

#    def sZx(self):

#b=Dynamics()
#b.Init(imagine=1)
#xi=b.rzero(force='real')
#print xi
#t,x,u=b.Traj(xi,U=1,Plot=1)
#b.nTraj(0,hold=1,U=1,imagine=1,RunList=[[xi[i]+1j*1./4096.*(random()-1/2) for i #in range(3)] for k in range(1)])
