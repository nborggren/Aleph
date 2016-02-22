"""The Algebra class serves as a backbone for basic operations, 
addition and multiplication, 
a torch to the abyss of differential equations ahead.
"""

__author__    = "Nathan A. Borggren (nborggren@gmail.com)"
__date__      = "2011-01-04"
__copyright__ = "Copyright (C) 2011 Nathan A. Borggren"
__license__   = "GNU LGPL Version 2.1"

from ROOT import TF3,TH1F,TFile, gDirectory
from numpy import *
from array import array
from swiginac import *  
from os import environ
from os import system as MySys
from time import time, strftime, gmtime
#from datetime import date
from Analysis import *

class Algebra:
    """ The Algebra Class 
    """
    ## Initialization
    def Init(self,dof=3,mon=5,S=False,system='Lorenz',
             f=["-sigma*x[0]+sigma*x[1]",
                "rho*x[0]-x[1]-x[0]*x[2]",
                "-beta*x[2]+x[0]*x[1]"],
             u=["1","y/x","x/y","x*z/y","x*y/z"],
             par={'sigma':10,'rho':28.,'beta':8./3.},
             spar = [symbol("\\sigma"),symbol("\\rho"),symbol("\\beta")],
             A=matrix([[-10.,10.,0,0,0],
                       [-1,0,28.,-1,0],
                       [-8./3.,0,0,0,1]]),
             B=matrix([[0,0,0],
                       [-1,1,0],
                       [1,-1,0],
                       [1,-1,1],
                       [1,1,-1]]),
             imagine=0,
             MyLog=1):

        self.dof = dof  
        self.mon = mon         
        self.S=S
        self.f=f
        self.u=u
        self.par=par
        self.spar=spar
        self.MyLog=MyLog

        if S==True:
            self.U = matrix([[symbol("u_{"+str(i)+"}")] for i in range(mon)])
            self.V = matrix([[symbol("v_{"+str(i)+"}")] for i in range(mon)])
            self.x = matrix([[symbol("x_{"+str(i)+"}")] for i in range(dof)])
            
            if system=='Lorenz':
                self.A=matrix([[-spar[0],spar[0],0,0,0],
                               [-1,0,spar[1],-1,0],
                               [-spar[2],0,0,0,1]])
                
                self.B=matrix([[0,0,0],
                               [-1,1,0],
                               [1,-1,0],
                               [1,-1,1],
                               [1,1,-1]])

                self.M = self.GetM()

        else:
            self.A = A
            self.B = B
            self.M = self.GetM()

        self.media = {'image':[],
                      'table':[],
                      'eqn':[],
                      'video':[],
                      'dot':[]}

        self.logbook = {'notes':[],
                        'inits':[]}

        self.root = [] 
        self.dtime= str(time()).split('.')[0] #date.isoformat(date.today())
        self.current,self.current2 = cname(0),cname(0)
        self.name = strftime("%Y_%b_%d.%H_%M", gmtime())+"_"
        self.output = self.name+"_"
        self.scurrent,self.scurrent2 = self.current.next(), self.current2.next()

    ## Calulates matrix product \f$M=BA\f$
    def GetM(self):
        return self.B*self.A

    ## Calulates structure constants in \f$u_i*u_j=\frac{1}{2}(\delta_{ik}M_{kj}+\delta_{ij}M_{jk})\f$ 
    def Struct(self,i,j,k):
        return (delta(i,k)*self.M[k,j]+delta(i,j)*self.M[j,k])/2

    ## the product of the ith and jth basis vector in the algebra \f$u*v\f$
    def mult(self,i,j):
        return [self.Struct(k,i,j) for k in range(self.mon)]

    ## the product of vectors in the algebra \f$u*v\f$
    def Mult(self,U,V):

        try:
            len(U)
        except TypeError:
            return [U*V[i] for i in range(self.mon)]

        try:
            tmp = [0 for i in range(self.mon)]
            for i in range(self.mon):
                for j in range(self.mon):
                    for k in range(self.mon):
                        tmp[k] = tmp[k]+U[i]*V[j]*self.Struct(k,i,j)
            return tmp

        except RuntimeError:
            tmp = [0 for i in range(self.mon)]
            for i in range(self.mon):
                for j in range(self.mon):
                    for k in range(self.mon):
                        tmp[k] = tmp[k]+U[i][0]*V[j][0]*self.Struct(k,i,j)
            return tmp
 
    ## the sum of vectors in the algebra \f$u+v\f$
    def Sum(self,U,V):
        try:
            return [U[i]+V[i] for i in range(self.mon)]

        except RuntimeError:
            return [U[i][0]+V[i][0] for i in range(self.mon)]

    ## we can call it additon too. \f$u+v\f$
    def Add(self,U,V):
        try:
            return [U[i]+V[i] for i in range(self.mon)]

        except RuntimeError:
            return [U[i][0]+V[i][0] for i in range(self.mon)]

    ## the associativity relations
    def Assoc(self,U,V,W):
        return [self.Mult(self.Mult(U,V),W)[i]-self.Mult(U,self.Mult(V,W))[i] for i in range(self.mon)]

    def IsJordan(self,U,V):
        return self.Mult(self.Mult(U,V),self.Mult(U,U))==self.Mult(U,self.Mult(V,self.Mult(U,U)))

    ## numeric or symbolic vector of point in monomial space
    def Mons(self,x):
        return [1,x[1]/x[0],x[0]/x[1],x[0]*x[2]/x[1],x[0]*x[1]/x[2]] 

    def tMons(self,x):
        tmp = [1,x[1]/x[0],x[0]/x[1],x[0]*x[2]/x[1],x[0]*x[1]/x[2]] 
        return [atan(j) for j in tmp]

    def PlotMons(self,U):
        sProp = ''
        for i in self.Sum(['['+str(i)+']*' for i in range(self.mon)],self.u):
            sProp=sProp+i+"+"
        return sProp[:-1]

    ## The time derivative of U comes from \f$\dot{U_i}=U_iM_{ij}U_j\f$,\f$\dot{x_i}= \sum_{j=1}^mA_{ij}\prod_{k=1}^nx_k^{B_{jk}}\f$
    def udot(self,U,x=False):
        if x==True:
            U = self.Mons(U)

        tmp = [sum([self.M[j,i]*U[i] for j in range(self.mon)]) for i in range(self.mon)]
     
        return [self.Mult(U[i],tmp[i]) for i in range(self.mon)]

    def Publish(self,name,title,path="./"):

        f=open(environ["ALEPHPATH"]+'lib/template.tex')
        g=open(name+'.tex','w')

        img_beg = ["\\begin{figure}[h]",
                   "\\begin{center}"]

        img_end = ["\\end{center}",
                   "\\end{figure}"]

        tbl_beg = ["\\begin{table}[htp]",
                   "\\centering",
                   "\\begin{tabular}{"+"*{"+str(self.mon+1)+"}{|>{\centering}p{"+str(15./(self.mon+1))[:3]+"cm}}|} ",
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

        for i in self.media['image']:
            for j in img_beg:
                IMAGES=IMAGES+j+" \n "

            IMAGES = IMAGES+ "\epsfig{file="+i+",scale = 0.8} \n"

            for j in img_end:
                IMAGES=IMAGES+j+" \n "

        for i in self.media['table']:
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

        for i in self.media['dot']:
            MySys("$DOTPATH./ladot "+i+".ladot")
            DOTS=DOTS+"\\input{"+i+".tex} \n"
            DOTS=DOTS+"\\includegraphics[width=6.5in,height=8in]{"+i+".ps} \n"

        for i in self.media['eqn']:
            for j in eqn_beg:
                EQNS=EQNS+j+" \n "

            EQNS = EQNS+i+" \n"

            for j in eqn_end:
                EQNS=EQNS+j+" \n "

        for i in self.media['video']:
            for j in vid_beg:
                VIDEOS=VIDEOS+j+" \n "

            VIDEOS = VIDEOS + "text={\\small("+i+")}]{10cm}{10cm}{"+i+"} \n "

            for j in vid_end:
                VIDEOS=VIDEOS+j+" \n "

        for line in f:
            line=line.replace('NOTES',title)
            line=line.replace('IMAGES',IMAGES)
            line=line.replace('TABLES',TABLES)
            line=line.replace('EQNS',EQNS)
            line=line.replace('DOTS',DOTS)
            line=line.replace('VIDEOS',VIDEOS)
            g.write(line)

        g.close()
                    
        MySys("latex "+name+".tex")
        MySys("dvipdf "+name+".dvi")

    ## returns TeX for the multiplication table
    def mtable(self):
        labels = ["$U_"+str(i+1)+"$" for i in range(self.mon)]
 
        output = "   "
        for i in range(self.mon):
            output = output+"$U_"+str(i+1)+"$ & $"
            for j in range(self.mon):
                tmp=self.mult(i,j)
                #print tmp
                for k in range(self.mon):
                    if tmp[k]!=0:
                        output = output + "("+str(tmp[k])+")U_"+str(k+1)+"+"
                if output[-1]=="+":
                    output=output[:-1]+ "$ &$"
                else:
                    output = output + " 0 $ &$"

            output = output[:-2] + "\n \\tabularnewline \n \hline \n"
                
        return output

    ##BasicGraph, default conneted graph of 
    def BoseHubbard(q,name='demo',k = 0, conds = [],outs = [],ps=1):
        
        points = [i for i in range(q)]

        f = open(name+'.dot','w')
        f.write('digraph ' + name + '{ \n')
        f.write('       edge [dir=none] \n')

        for j in points:
            if i<j:                   #+1 line for double counting
                f.write(str(i)+'->'+str(j)+'; \n ')
        f.write('}')
        f.close()

        MySys("dot -Tps "+name+".dot > "+name+".ps")

        return
    ## makes a ladot file demonstrating the recursion relations 
    def Ladot(self,name,D=False):
        f = open(name+".ladot",'w')
        ltr = [chr(i) for i in xrange(ord('i'), ord('i')+self.dof)]
        alp = [chr(i) for i in xrange(ord('a'), ord('a')+self.mon)] 
        f.write("graph "+name+" { \n \n    node [shape=box,width=2.5, color = \"blue\"] \n")
        sym = [symbol(i) for i in ltr]
        a=matrix([[zz+1 for zz in sym]])
        b=a*self.A

        for ii in range(self.mon):
            tmp = " "*4+alp[ii]+" [label = \"$p_{"
            tmp2 = "["+str(b[0,ii])+"]$(6)\"]"
            for jj in range(self.dof):
                if self.B[ii,jj]==0:
                    tmp = tmp+ltr[jj]
                elif self.B[ii,jj]>0:
                    tmp = tmp+ltr[jj]+str(-self.B[ii,jj])
                else:
                    tmp = tmp+ltr[jj]+"+"+str(-self.B[ii,jj])

            tmp = tmp+"}"
            f.write(tmp+tmp2+" \n")

        f.write("    edge [color = \"red\"] \n \n")

        for ii in alp:
            for jj in alp:
                if ii<jj:
                    f.write("    "+ii+" -- "+jj+"; \n")

        f.write("} \n")
        self.dot.append(name)

    ##Return a root TF2 object for plotting vector u
    def GetTF2(self,u):
        if self.dof!=2: 
            return "sorry"
        fnc = ""
        for i in range(self.mon-1):
            fnc=fnc+str(u[i])+"*"+self.U[i]+"+"
        fnc=fnc+str(u[self.mon-1])+"*"+self.U[self.mon-1]
        return TF2("sys",fnc,0,10,0,10)

    ##Return a root TF2 object plotting vector (over monomials) u 
    def GetTF3(self,u):
        if self.dof!=3: 
            return "sorry"
        fnc = ""
        for i in range(self.mon-1):
            fnc=fnc+str(u[i])+"*"+self.u[i]+"+"
        fnc=fnc+str(u[self.mon-1])+"*"+self.u[self.mon-1]
        return TF3("lorenz",fnc,-30,30,-30,30,0,50)

    ##Newton solver to find zeroes of F, (Null Space of M)
    def GetZeroes(self):
        return 

    def write(self,output):
        zwom = TFile(output+".root","recreate")
        for i in self.root:
            try:
                i.Write()
            except AttributeError:
                continue

        #
        #for i in self.hU:i.Write()
        #for i in self.tU:i.Write()
        #for i in self.hSect:i.Write()

## delta function returns true if both inputs are the same.
def delta(i,j):
    return float(i==j)

def MyPow(Dyn,u,m):
    if m == 1:
        return u
    else:
        return Dyn.Mult(u,MyPow(Dyn,u,m-1))

##The Algebra system class.
#A quasimonomial dynamical system can be written for i from 1 to n (dof) as 
#\f$\dot{x_i}=x_i\sum_{j=1}^mA_{ij}\prod_{k=1}^nx_k^{B_{jk}}\f$ \cite{Figuerdo} 


            

