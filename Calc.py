from Tkinter import *
from ROOT import *

class Calc(Frame):

    def Start(self, master=None):
        Frame.__init__(self,master)
        self.grid(sticky='nesw')

        self.sOps = ['traj','ntraj','hist','lang','nlang','mesh','meshtraj','myFFT','Reset']

        self.dOps = ['run','fpt']

        self.Widgets(label=False)

    ##Widgets for Gui
    def Widgets(self, label=False):

        top=self.winfo_toplevel()                
        top.rowconfigure(0, weight=1)            
        top.columnconfigure(0, weight=1)         
        self.rowconfigure(0, weight=1)           
        self.columnconfigure(0, weight=1)

        #?????????????????????????????????????????????????????????????????
        #Algebra widgets
        #?????????????????????????????????????????????????????????????????
 
        #frame displaying velocity field
        self.wForce = Frame(self)
        self.wForce.grid(row=2,rowspan=4,column=1,columnspan=8,sticky='nesw')
        self.wF = []
        for i in range(self.dof):
            print self.f[i]
            if label==True:
                self.wF.append(Label(self.wForce,text="f_"+str(i)+"= "+self.f[i]))
                self.wF[-1].grid(row=2*i,rowspan=2,columnspan=8)
            else:
                self.wF.append(Label(self.wForce,text="f_"+str(i)))
                self.wF[-1].grid(row=2*i,rowspan=2,columnspan=8)

        #frame displaying parameters of velocity field
        self.wPar = Frame(self)
        self.wPar.grid(row=2, column=10, columnspan=8,rowspan=4,sticky='nesw')

        self.sPar = [StringVar() for i in range(len(self.par))]
        for i in range(len(self.par)):
            tmp = self.par.keys()[i] 
            self.sPar[i].set(str(self.par[tmp]))
        self.parInits = []
        self.lInits = []

        nn = len(self.par.values())
        for i in range(nn):
            self.parInits.append(Entry(self.wPar,textvariable=self.sPar[i]))
            self.lInits.append(Label(self.wPar,text=self.par.keys()[i]))
            self.parInits[-1].grid(row =i/5+1, column=6+2*(i%5)+1)
            self.lInits[-1].grid(row =i/5+1, column=6+2*(i%5))

        #frame displaying parameters for trajectories and langevin trajectories

        self.tpar = {"t0":0.,
                     "t1":25000.,
                     "t2":5000,
                     "t3":1,
                     "t4":50}

        self.wtPar = Frame(self)
        self.wtPar.grid(row=25, column=10, columnspan=8,rowspan=4,sticky='nesw')

        self.stPar = [StringVar() for i in range(len(self.tpar))]
        for i in range(len(self.tpar)):
            tmp = self.tpar.keys()[i] 
            self.stPar[i].set(str(self.tpar[tmp]))
        self.tparInits = []
        self.tlInits = []

        for i in range(len(self.tpar.values())):
            self.tparInits.append(Entry(self.wtPar,textvariable=self.stPar[i]))
            self.tparInits[-1].grid(row = 1, column =6+2*(i%5)+1)
            self.tlInits.append(Label(self.wtPar,text=self.tpar.keys()[i]))
            self.tlInits[-1].grid(row = 1, column=6+2*(i%5))

        #frame displaying parameters of diffusion tensor

        self.wdPar = Frame(self)
        self.wdPar.grid(row=15, column=10, columnspan=8,rowspan=2,sticky='nesw')


        self.sdPar = [StringVar() for i in range(len(self.dpar))]
        for i in range(len(self.dpar)):
            tmp = self.dpar.keys()[i] 
            self.sdPar[i].set(str(self.dpar[tmp]))
        self.dparInits = []
        self.dlInits = []

        for i in range(len(self.dpar.values())):
            self.dparInits.append(Entry(self.wdPar,textvariable=self.sdPar[i]))
            self.dparInits[-1].grid(row = 1, column =6+2*(i%5)+1)
            self.dlInits.append(Label(self.wdPar,text=self.dpar.keys()[i]))
            self.dlInits[-1].grid(row = 1, column=6+2*(i%5))

        #frame displaying noise terms for langevin equations

        self.lngPar = Frame(self)
        self.lngPar.grid(row=17, column=10, columnspan=8,rowspan=2,sticky='nesw')
        self.slngPar = [StringVar() for i in range(self.dof)]
        self.slngInits = []

        for i in range(self.dof):
            self.slngPar[i].set(self.nexpr[i])
        for i in range(self.dof):
            self.slngInits.append(Entry(self.lngPar,textvariable=self.slngPar[i]))
            self.slngInits[i].grid(row = 1, column =6+i)

        ##collect images used for buttons
        self.phX = [PhotoImage(file="/home/nborggren/Aleph/lib/micros/x"+str(i)+".gif") for i in range(self.dof)]
        self.wInits = Frame(self)
        self.wInits.grid(row=8,rowspan=4,column=1,columnspan=2,sticky='nesw')
        self.xzero = [StringVar() for i in range(self.dof)]
        self.suzero = [StringVar() for i in range(self.mon)]

        tmp = self.rzero()
        utmp = self.Mons(tmp)
        sx = ['x_'+str(i) for i in range(self.dof)]
        for i in range(self.mon):
            self.suzero[i].set(str(utmp[i])[:6])
        for i in range(self.dof):
            self.xzero[i].set(str(tmp[i]))
        self.sxInits = [Entry(self.wInits,textvariable=self.xzero[i]) for i in range(self.dof)]
        self.bxInits = [Button(self.wInits,image=self.phX[i]) for i in range(self.dof)]
        ##scale for initial conditions
        for i in range(self.dof):
            self.sxInits[i].grid(row=i,column=1) 
            exec("self.bxInits["+str(i)+"].bind(\"<Button-1>\", self.ev_x"+str(i)+')')
        ##buttons for form generators
            self.bxInits[i].grid(row=i,column=2) 

        #soul operators
        self.wOps = Frame(self)
        self.wOps.grid(row=8, column=14, columnspan=3,rowspan=6,sticky='nesw')
        self.aOps = []
        for i in range(len(self.sOps)):
            print i
            exec("self.aOps.append(Button(self.wOps,text=\""+self.sOps[i]+"\",command=self.ev_"+self.sOps[i]+"))")
            self.aOps[i].grid(row=i/3,column=i%3)

        #soul operators
        self.wdOps = Frame(self)
        self.wdOps.grid(row=19, column=14, columnspan=5,rowspan=6,sticky='nesw')
        self.adOps = []
        for i in range(len(self.dOps)):
            print i
            exec("self.adOps.append(Button(self.wdOps,text=\""+self.dOps[i]+"\",command=self.ev_"+self.dOps[i]+"))")
            self.adOps[i].grid(row=i/3,column=i%3)

        lu0 = Label(self.wdOps,text='u0')
        lu0.grid(row=2,column=0)
        self.su0 = StringVar()
        self.su0.set("1./pow(50*2.5066,2)*exp((-pow(x[0]-"+self.xzero[0].get()+",2)-pow(x[1]-"+self.xzero[1].get()+",2))/(2*pow(50,2)))")
        wu0 = Entry(self.wdOps,textvariable=self.su0)
        wu0.grid(row=2,column=1)

        lname = Label(self.wdOps,text='name')
        lname.grid(row=3,column=0)
        self.soutput = StringVar()
        self.soutput.set(self.output)
        wname = Entry(self.wdOps,textvariable=self.soutput)
        wname.grid(row=3,column=1)
 
        zouts = ['lies','maybe','truth']
        self.wExpr = Frame(self)
        self.wExpr.grid(row=16, column=1, columnspan=8,rowspan=3,sticky='nesw')
        self.var = [Label(self.wExpr,text=zouts[i]) for i in range(3)]
        self.remarks = [StringVar() for i in range(3)]
        self.aExpr = [Label(self.wExpr,textvariable=self.remarks[i]) for i in range(3)]
        for i in range(3):
            self.aExpr[i].grid(row=i,column=1)
            self.var[i].grid(row=i,column=10)

        self.phPar = [PhotoImage(file="/home/nborggren/Aleph/lib/micros/"+i+".gif") for i in ['sigma','rho','beta']]
        self.wPar = Frame(self)
        #self.wPar.grid(row=12,column=1,columnspan = 3)
        self.bPar = [Button(self.wPar,image=i) for i in self.phPar]
        #for i in range(3):self.bPar[i].grid(row=0,column=i)
         
    def ev_Reset(self):

        for i in range(len(self.sPar)):
            print self.par.keys()[i]
            self.par[self.par.keys()[i]]=float(self.sPar[i].get())

        #no se porque
        for i in range(len(self.sPar)):
            print self.par.keys()[i]
            self.par[self.par.keys()[i]]=float(self.sPar[i].get())

        for i in range(len(self.sdPar)):
            self.dpar[self.dpar.keys()[i]]=float(self.sdPar[i].get())
        
        try:
            x = [float(self.xzero[i].get()) for i in range(self.dof)]

        except ValueError:
            x = [complex(self.xzero[i].get()) for i in range(self.dof)]
            self.imagine==1
            #self.Traj(x,U=1,Plot=1)
            #return
        
        self.reset()
       
    def current(self,ohwell):
        tmp = [self.xzero[i].get() for i in range(self.dof)]
        utmp = self.Mons(tmp)       
        for i in range(self.mon):
            self.suzero[i] = str(utmp[i])[:6]
       
        pass

    def ev_traj(self):
        try:
            x = [float(self.xzero[i].get()) for i in range(self.dof)]

        except ValueError:
            x = [complex(self.xzero[i].get()) for i in range(self.dof)]
            self.imagine==1

        self.Traj(x,U=1,Plot=1)

        return

    def ev_ntraj(self,Plot=1):
        self.nTraj(50,U=1)
        pass

    def ev_lang(self):
        x = [float(self.xzero[i].get()) for i in range(self.dof)]
        print 'langevining'
        
        for i in range(len(self.sPar)):
            
            self.par[self.par.keys()[i]]=float(self.sPar[i].get())
            
        self.reset()
        try:
            #self.Lang(x,m=2,wt=float(self.sdPar[2].get()),Plot=1)
            self.Lang(x,m=2,wt=self.par['Dxx'],Plot=1) #buggy fixed weight
        except ValueError:
            self.Lang(x,m=2,wt=100.,Plot=1)
        pass

    def ev_nlang(self,Plot=1):
        for i in range(self.dof):
            self.par[self.par.keys()[i]]=float(self.sPar[i].get())
        self.reset()
        self.nLang(50,Plot=1)
        pass

    def ev_hist(self):
        c3 = TCanvas()
        self.hist.Draw("colz")
        zwom=input("sheat")
        pass

    def ev_mesh(self):
        self.Mesh(self.output+".xml")
        pass

    def ev_meshtraj(self):
        #for i in range(self.dof):
        #    print self.sPar[i].get()
        #    self.par[self.par.keys()[i]]=self.sPar[i].get()
        #self.ev_Reset()
        self.MeshTraj(int(self.tsteps),1,self.dt,sz=15,boundary=1)
        pass

    def ev_run(self):

        #self.output=self.soutput.get()
        #print self.output
        self.Run(self.output)

    def ev_fpt(self):
        self.FirstPassage()
        pass




