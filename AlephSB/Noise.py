#???????????????????????????????????????????????????????????????????????????????????
##@package Noise
#Classes for Dynamical Systems with random force from ROOT histograms. 
#@author Nathan Borggren
#
#???????????????????????????????????????????????????????????????????????????????????

from Dynamics import *
from random import choice
from cmath import sqrt as iSqrt

##Noise class studies the dynamics of a system with a random force (default gaussian) 
class Noise(Dynamics):
    def SetNoise(self,Expr="gaus(0)",Fmin=-75,Fmax=75,wt=1.,imagine=1):
        self.wt = wt
        Dynamics.Init(self,imagine=imagine)
        self.noise = [TF1("noise",Expr,Fmin,Fmax) for i in range(self.dof)]
        self.hnoise = [TH1F("hnoise_"+str(i),"hnoise_"+str(i),100,Fmin,Fmax) for i in range(self.dof)]
        for i in self.noise:i.SetParameters(300.,0.,wt)

    ##Finite difference integration step calling a random number from the noise expression. 
    def NoiseStep(self,x,dt):

        wt = self.wt

        if self.imagine==1:
            #print x
            tmp = [self.noise[i].GetRandom()+self.noise[i].GetRandom()*iSqrt(-1) for i in range(self.dof)]
            for i in range(self.dof):
                self.hnoise[i].Fill(wt*tmp[i].real)
        else:
            tmp = [self.noise[i].GetRandom() for i in range(self.dof)]

        #for i in range(self.dof):
        #    self.hnoise[i].Fill(wt*tmp[i])

	xr = [x[i]+wt*tmp[i] for i in range(self.dof)]
	a1 = self.Flow(xr)
 
        return [x[i] + dt*a1[i] for i in range(self.dof)] 
 
    ##Returns a trajectory of langevin dynamics from input initial conditions.
    def Lang(self,x,m=50,Fill=None,U=None,isPos=0,wt=1.,Plot=0,fast=1):

        traj=[array('d') for i in range(self.dof)]
        xraj=array('d')   
        xraj.append(self.ti)

        if Fill==1:
            ost = ""
            for j in range(self.dof):
                ost = ost+"x["+str(j)+"],"
            ost = ost[:-1]

        if U==1:
            utraj=[array('d') for i in range(self.mon)]
            tmp=self.Mons(x)
            for i in range(self.mon):
                utraj[i].append(tmp[i])

        for i in range(self.dof):traj[i].append(x[i])
       
        for j in range(self.tsteps):
            x = self.Noisestep(x,self.dt,wt=wt)
            xraj.append((j+1)*self.dt+self.ti)
            for i in range(self.dof):traj[i].append(x[i])
            
            if U==1:
                tmp=self.Mons(x)
                for i in range(self.mon):utraj[i].append(tmp[i])

            if Fill==1 and j>self.roc*self.tsteps:
                exec("self.hist.Fill("+ost+")")
                for i in range(len(self.pSect)):
                    exec(self.cSect[i])
                    if cond:
                        exec(self.fSect[i])  #must define fill array
                        exec("self.hSect[i].Fill("+str(fill[0])+','+str(fill[1])+")")
        
        self.last=x

        print len(xraj),len(traj),"here we are"
        if Plot==1:
            if U==1:
                if self.imagine==1:
                    c1,c2,c4,c5=self.Prettify(xraj,traj,uu=utraj)
                else:
                    c1,c2,c4=self.Prettify(xraj,traj,uu=utraj)
            else:
                c2,c4=self.Prettify(xraj,traj)

            c2.Draw()
            c4.Draw()
            zwom = input("sheat")
        if U ==1:
            return xraj, traj, utraj

        else:
            return xraj, traj

        ##Does n such trajectories.  Good for initial condition guess for diffusion equations.
    def nLang(self,n,m=50,wt=1.):
        try:
            self.hist
        except AttributeError:
            self.hReset()
        system("echo "+str(n)+" trajectories started at: " +strftime('%y.%m.%d -- %H.%M.%S ')+" >> "+self.output+".log")
        t=0
        before=time()
        
        for i in range(n):
            if i%50==0:
                t+=(time()-before)/60.
                before = time()
                print t, " minutes for ", i, "trajectories"
            x=self.rzero()
            self.Lang(x,m,Fill=1,wt=wt)

        system("echo "+str(n)+" trajectories finished at: " +strftime('%y.%m.%d -- %H.%M.%S ')+" >> "+self.output+".log")
        return
