import os, sys
sys.path.append(os.path.join('/usr/local/bin/root/','lib'))

from Tkinter import *
from ROOT import *

class Dissonance:
    def __init__(self, root):



        self.myParent = root

        self.top_frame = Frame(root)
        self.top_frame.pack(side=RIGHT,fill=BOTH,expand=YES)

        self.left_frame = Frame(self.top_frame, background="grey",
                                borderwidth=5,  relief=RIDGE,
                                height=350, 
                                width=650, 
                                ) ###
        self.left_frame.pack(side=RIGHT,
                             fill=BOTH, 
                             expand=YES,
                             )  ###
        
        self.micros = Frame(root)
        self.micros.pack()

        self.rght = Frame(root)
        self.rght.pack()

        self.bJos = Button(self.micros)
        self.bRes = Button(self.micros)

        photo = PhotoImage(file="./lib/micros/jos.GIF")
        self.bJos = Label(self.micros, image=photo)
        self.bJos.photo = photo
        self.bJos["background"] = "blue"
        self.bJos.bind("<Button-1>", self.bJosClick) 
        self.bJos.pack()

        photo1 = PhotoImage(file="./lib/micros/res.GIF")
        self.bRes = Label(self.micros, image=photo1)
        self.bRes.photo = photo1
        self.bRes["background"] = "blue"
        self.bRes.bind("<Button-1>", self.bResClick) 
        self.bRes.pack()

        
        photo2 = PhotoImage(file="./lib/micros/inductor.GIF")
        self.bInd = Label(self.micros, image=photo2)
        self.bInd.photo = photo2
        self.bInd["background"] = "blue"
        self.bInd.bind("<Button-1>", self.bIndClick) 
        self.bInd.pack()

        photo3 = PhotoImage(file="./lib/micros/cap.GIF")
        self.bCap = Label(self.micros, image=photo3)
        self.bCap.photo = photo3
        self.bCap["background"] = "blue"
        self.bCap.bind("<Button-1>", self.bCapClick) 
        self.bCap.pack()

        photo4 = PhotoImage(file="./lib/micros/source.GIF")
        self.bSrc = Label(self.micros, image=photo4)
        self.bSrc.photo = photo4
        self.bSrc["background"] = "blue"
        self.bSrc.bind("<Button-1>", self.bSrcClick) 
        self.bSrc.pack()

        photo5 = PhotoImage(file="./lib/micros/sink.GIF")
        self.bSnk = Label(self.micros, image=photo5)
        self.bSnk.photo = photo5
        self.bSnk["background"] = "blue"
        self.bSnk.bind("<Button-1>", self.bSnkClick) 
        self.bSnk.pack()

        photo6 = PhotoImage(file="./lib/micros/gnd.GIF")
        self.bGnd = Label(self.micros, image=photo6)
        self.bGnd.photo = photo6
        self.bGnd["background"] = "blue"
        self.bGnd.bind("<Button-1>", self.bGndClick) 
        self.bGnd.pack()

        self.Img = Canvas(self.left_frame,width=1000,height=500)
        self.Img.pack()
        self.Img.create_line(0,100,470,0,fill="red",dash=(4,4))

    def bJosClick(self, event):
        if self.bJos["background"] in ["blue","red"]: 
            self.bJos["background"] = "yellow"
        elif self.bJos["background"] == "yellow":
            self.bJos["background"] = "green"

        else:
            self.bJos["background"] = "red"

    def bResClick(self, event):

        if self.bRes["background"] in ["blue","red"]: 
            self.bRes["background"] = "yellow"

        elif self.bRes["background"] == "yellow":
            self.bRes["background"] = "green"
        else:
            self.bRes["background"] = "red"

    def bIndClick(self, event):

        if self.bInd["background"] in ["blue","red"]: 
            self.bInd["background"] = "yellow"

        elif self.bInd["background"] == "yellow":
            self.bInd["background"] = "green"
        else:
            self.bInd["background"] = "red"

    def bCapClick(self, event):

        if self.bCap["background"] in ["blue","red"]: 
            self.bCap["background"] = "yellow"

        elif self.bCap["background"] == "yellow":
            self.bCap["background"] = "green"
        else:
            self.bCap["background"] = "red"

    def bSrcClick(self, event):

        if self.bSrc["background"] in ["blue","red"]: 
            self.bSrc["background"] = "yellow"

        elif self.bSrc["background"] == "yellow":
            self.bSrc["background"] = "green"
        else:
            self.bSrc["background"] = "red"

    def bSnkClick(self, event):

        if self.bSnk["background"] in ["blue","red"]: 
            self.bSnk["background"] = "yellow"

        elif self.bSnk["background"] == "yellow":
            self.bSnk["background"] = "green"
        else:
            self.bSnk["background"] = "red"  

    def bGndClick(self, event):

        if self.bGnd["background"] in ["blue","red"]: 
            self.bGnd["background"] = "yellow"

        elif self.bGnd["background"] == "yellow":
            self.bGnd["background"] = "green"
        else:
            self.bGnd["background"] = "red"  


class Jos:

    IxRm = 26.0       #);   Subgap voltage Ic*Rm in mV.
    Jc   = 4500.0     #); Critical current density in A/cm2.
    RmRn = 20.0       #);   Dimensionless ratio of subgap (Rm) and normal (Rn) resistances.

    #I = Is(phi)+In(V)+Id(Vdot)+If(t) 
    # pass
    
#?????????????????????????????????????????????????????????????????????????????????????????????????????
#?Unit Conversion to Pscan
#?
#?
#?????????????????????????????????????????????????????????????????????????????????????????????????????

def Convert_to_Pscan(type,num):
    if type == 'R':
        return num/2.38         #One PSCAN resistance unit corresponds to 2.38 Ohms.
    elif type == 'L':
        return num/2.64         #One PSCAN unit for inductance corresponds to 2.64 pH.
    elif type == 'I':
        return num/0.125        #One PSCAN unit for current corresponds to 0.125 mA.
    #elif type == 'V':
    #    return num/2.6          #For calculation of bias current it is supposed that biasVolt=2.6 mV
    #query http://pavel.physics.sunysb.edu/RSFQ/Lib/units.html
    else:
        return 'enter R for resistance, L for inductance, I for Current.'

def Convert_away_Pscan(type,num):
    if type == 'R':
        return num*2.38         #One PSCAN resistance unit corresponds to 2.38 Ohms.
    elif type == 'L':
        return num*2.64         #One PSCAN unit for inductance corresponds to 2.64 pH.
    elif type == 'I':
        return num*0.125        #One PSCAN unit for current corresponds to 0.125 mA.
    #query http://pavel.physics.sunysb.edu/RSFQ/Lib/units.html
    else:
        return 'enter R for resistance, L for inductance, I for Current.'
  
    


root = Tk()

dis = Dissonance(root)
root.mainloop()

