from scanner.labin import specscope,webocean,uniocean
from scanner import labrc as rc

class specscope3d(specscope):
        scale=1
        gxmax=200
        gymax=180
        gzmax=30
        gzmin=4

        def awrite(self,comm):
            self.ard.write((comm+"\r\n").encode())
            
        def __init__(self, *args, **kwargs):
                specscope.__init__(self, *args, **kwargs)

        def adrive(self,peri=10,gap=0):
                import serial
                if not hasattr(self,'ard'):
                    try:
                        self.ard=serial.Serial(rc.ard_port,baudrate=rc.ard_speed,timeout=2)
                    except:
                        print("cannot connect to Velleman")
                        return -1
                self.awrite('G21') #milimetres
                self.awrite('M17') #START motors
                self.awrite('G90')
                return self.acomm()

        def acomm(self,isok=False):
                rep=[self.ard.readline()]
                while (len(rep[-1])>0):
                        rep.append(self.ard.readline())
                        if isok and rep[-1][:2]=='ok': break
                return [r for r in rep[:-1] if len(r.strip())>0]

        def ahome(self):
                self.awrite("G28 X")
                self.awrite("G28 Y")
                if hasattr(rc,'xy_cent'):
                    self.goto(*rc.xy_cent)
                return self.acomm()

        def rate(self,ch=1,nstep=100,wait=0,ain=0,ntries=3,loud=0):
                if ch==1:
                        self.gx+=nstep*self.scale
                        self.awrite("G1 X%.1f"%self.gx)
                elif ch==2:
                        self.gy+=nstep*self.scale
                        self.awrite("G1 Y%.1f"%self.gy)
                return self.acomm(True)

        def goto(self,xpos,ypos):
                if xpos>=0 and xpos<=self.gxmax: self.gx=xpos
                if ypos>=0 and ypos<=self.gymax: self.gy=ypos
                self.awrite("G1 X%.1f Y%.1f"%(self.gx,self.gy))
                return self.acomm(True)

        def focus(self,zpos=10):
            
            prof=[]
            self.setup(prange=[200,400],integ=10,aver=10)
            while (zpos>self.gzmin) && (zpos<self.gzmax):
                self.awrite("G1 Z%.1f"%zpos)
                self.measure()
                prof.append([zpos,self.last.sum()])
                if prof[-1][1]<prof[0][1]*0.8: break
            import numpy as np
            prof=np.array(prof)
            idx=np.polyfit(prof[:,0],prof[:,1],4)
            ip=prof[:,1].argmax()
            fx=np.r_[prof[0][i-2]:prof[0][i+3]:100j]
            self.zbest=fx[np.polyval(idx,fx).argmax()]
            return self.acom()

        def round(self,step=10,nstep=20,cent=[95,125],rad=75,run=0):
            #  from scanner.camera_anal import circum # stred ze 3 bodu
            import numpy as np
            ix,iy=np.indices((nstep,nstep))
            cx,cy=cent
            iy[::2]=iy[::2,::-1]
            sel=(ix*step-cx)**2+(iy*step-cy)**2<rad**2
            pos=np.r_[[ix[sel]*step,iy[sel]*step]].T
            self.omap={}
            if run>0:
                for p in pos:
                    self.goto(p[0],p[1])
                    pident="p%03i-%03i"%(p[0],p[1])
                    if run>1: self.omap[pident]=self.result()
            else:
                return pos
                
            def end(self):
                if hasattr(self,"ard") and self.ard!=None: 
                    self.ard.close()
                    del self.ard
                    
from scanner.comband import gettrans
#from scanner.labin import webfetch

class webocean3d(uniocean,specscope3d):

        def __init__(self, *args, **kwargs):
            if rc.java_ooijar.find("Serial")<0:
                rc.java_ooijar+=";"+rc.java_ooipath+"jSerialComm.jar"
            uniocean.__init__(self, *args, **kwargs)
            rep=self.adrive()
            #print(rep[-1])

        def fastline(self,nstep,dstep=10,axe='X',mchan=0,bin=1,imin=None,imax=None):
            import numpy as np
            query=self.makequery(chan=mchan)+"&step={}&dstep={}{}".format(nstep,axe,dstep)
            if bin>1: query+="&bin={}".format(bin)
            print("running "+query)
            self.ard.close()
            out=self.getresp(query)
            self.ard=serial.Serial(rc.ard_port,baudrate=rc.ard_speed,timeout=2)
            pos=np.array([int(c.text) for c in out.findall("pos")])
            data=np.array([[float(a) for a in c.text.split()] for c in out.findall("data")])
            avg=data[:,imin:imax].mean(1)
            return pos,avg
            
        def end(self):
            uniocean.end(self)
            specscope3d.end(self)
            