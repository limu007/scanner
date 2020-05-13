from .labin import specscope,webocean,uniocean,oceanjaz
from . import labrc as rc

class specscope3d(specscope):
        scale=1
        gxmax=200
        gymax=180
        gzmax=6
        gzmin=2.5

        def awrite(self,comm):
            self.ard.write((comm+"\r\n").encode())
            
        def __init__(self, *args, **kwargs):
                specscope.__init__(self, *args, **kwargs)
                self.gx=0
                self.gy=0
                self.gz=self.gzmin

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

        def ahome(self,pos=[]):
                self.awrite("G28 X")
                self.awrite("G28 Y")
                if len(pos)>1:
                    self.goto(*pos)
                elif hasattr(rc,'xy_cent'):
                    self.goto(*rc.xy_cent)
                return self.acomm()

        def rate(self,ch=1,nstep=100,wait=0,ain=0,ntries=3,loud=0):
                if ch==1:
                        self.gx+=nstep*self.scale
                        self.awrite("G1 X%.1f"%self.gx)
                elif ch==2:
                        self.gy+=nstep*self.scale
                        self.awrite("G1 Y%.1f"%self.gy)
                if loud>=0: return self.acomm(True)
                return []

        def goto(self,xpos,ypos,zpos=0):
                if self.gz<self.gzmin:
                    self.awrite("G1 X%.1f Y%.1f Z%.2f"%(self.gx,self.gy,self.gzmin))
                if xpos>=0 and xpos<=self.gxmax: self.gx=xpos
                if ypos>=0 and ypos<=self.gymax: self.gy=ypos
                if len(rc.xy_zslope)>2:
                    if zpos>0: 
                        self.gz=zpos
                        self.awrite("G1 X%.1f Y%.1f"%(self.gx,self.gy))
                        self.awrite("G1 X%.1f Y%.1f Z%.2f"%(self.gx,self.gy,self.gz))
                    else:
                        self.gz=rc.xy_zslope[0]+rc.xy_zslope[1]*self.gx+rc.xy_zslope[2]*self.gy
                        if len(rc.xy_zslope)>5:
                            self.gz+=rc.xy_zslope[3]*self.gx**2+rc.xy_zslope[4]*self.gy**2+rc.xy_zslope[5]*self.gx*self.gy
                        if self.gz<self.gzmin: self.gz=self.gzmin
                        if self.gz>self.gzmax: self.gz=self.gzmax
                        self.awrite("G1 X%.1f Y%.1f Z%.2f"%(self.gx,self.gy,self.gz))

                else:
                    self.awrite("G1 X%.2f Y%.2f"%(self.gx,self.gy))
                return self.acomm(True)

        def focus(self,zpos=0,zstep=.2,rep=0,integ=10,factval=[1,0,0],minval=1e5,bord=10):
            '''rep=1: sum all channels
               rep=2: individual channels
            '''
            prof=[]
            self.intime=integ
            import numpy as np
            #self.setup(prange=[1.2,2.2],integ=integ,aver=10)
            if len(factval)>0: self.intfact=factval
            while (zpos>self.gzmin) and (zpos<self.gzmax):
                self.awrite("G1 Z%.1f"%zpos)
                last=self.measure().copy()
                if np.iterable(self.dark): last-=self.dark
                if bord>0: last=last[:,bord:-bord]
                if last.sum()<minval: #no signal here
                    break
                if rep>1:
                    prof.append([zpos]+list(last.sum(1)))
                else:
                    prof.append([zpos,last.sum()])
                zpos+=zstep
                #if prof[-1][1]<prof[0][1]*0.8: break
            self.gz=self.gzmax
            prof=np.array(prof).T
            if rep>0 or len(prof[0])<5: return prof #no fitting
            k=1
            idx=np.polyfit(prof[0],prof[k],4)
            ip=prof[k].argmax()
            fx=np.r_[prof[0][ip-2]:prof[0][ip+3]:100j]
            self.zbest=fx[np.polyval(idx,fx).argmax()]
            return self.acom()

        def prepare(self):
            self.setup([1.1,2.3],10)
            self.dark=self.measure().copy() 
            self.ahome()
            #[np.roots(np.polyder(np.polyfit(prof[0],p,4))) for p in prof[1:3]]

        def ruler(self,axis=1,stp=1,nstp=60,ichan=0):
            import numpy as np
            side1=np.array([len(self.rate(axis,stp,loud=-1))+self.measure()[ichan].mean() for i in range(nstp)])
            self.rate(axis,-nstp*stp)
            splt=np.percentile(side1,[10,90]).mean()
            if splt/side1.max()>0.8: len1=nstp*stp
            else: len1=sum(side1>splt)*stp
            side2=np.array([len(self.rate(axis,-stp,loud=-1))+self.measure()[ichan].mean() for i in range(nstp)])
            self.rate(axis,nstp*stp)
            splt2=np.percentile(side2,[10,90]).mean()
            if splt2/side2.max()>0.8: len2=nstp*stp # mozna neco lepsiho?
            else: len2=sum(side2>splt2)*stp
            return len1,len2

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

def zcalib(self,pos,factval=[1,1,0]):
    import numpy as np
    self.goto(pos[0],pos[1])
    prof=self.focus(self.gzmin+0.1,0.1,rep=2,integ=10,factval=factval,minval=10,bord=100)
    idx=[np.polyfit(prof[0],p,4) for p in prof[1:sum(factval)+1]]
    roughpos=[prof[0][p.argmax()] for p in prof[1:sum(factval)+1]]
    pox=[np.roots(np.polyder(p)) for p in idx]
    vax=[np.polyval(idx[i],pox[i]).real.argmax() for i in range(len(pox))]
    maxpos=[pox[i][vax[i]].real for i in range(len(pox))]
    chis=[sum((prof[i+1]-np.polyval(idx[i],prof[0]))**2)/prof[i+1].var()/len(prof[0]) for i in range(len(idx))]
    sigs=pox=[1/np.polyval(np.polyder(np.polyder(idx[i])),maxpos[i]) for i in range(len(idx))]
    return maxpos,chis,sigs,roughpos

from scanner.comband import gettrans
#from scanner.labin import webfetch

class ocean3d(oceanjaz,specscope3d):

        def __init__(self, *args, **kwargs):
            #if rc.java_ooijar.find("Serial")<0:
            #    rc.java_ooijar+=";"+rc.java_ooipath+"jSerialComm.jar"
            oceanjaz.__init__(self, *args, **kwargs)
            if self.dsize<=0: return None # init failed
            #if len(args)>1: return 
            rep=self.adrive()
            #print(rep[-1])
            self.gx=0
            self.gy=0
            self.gz=self.gzmin


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
            if hasattr(self,"ard") and self.ard!=None: 
                self.ard.close()
                del self.ard
            oceanjaz.end(self)
            #specscope3d.end(self)

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


## historical scanning procedure
def scan_save(self,xstep=-500,ystep=500,xcnt=50,ycnt=50,radius=0,oname='C:\\Users\\Lenovo\\Documents\\Data\\quad_line',
            center=False,format='txt',control=None,swapaxes=False,centerfirst=False,debug=False,returnlast=False,recalib=0):
        ''' 2-d scanning of the sample - saving each row in separate file (given by oname)
        if radius>0: X-Y coverage of circular wafer
            now radius in stepper counts
        if swapaxes: scanning second dir. first
        if centerfirst: just one point in first line (center)
        if returnlast: return to original position after scan
        '''
        from numpy import array,save,zeros,savetxt
        from math import sqrt
        import sys
        from time import sleep

        #global gx,gy
        self.gx,self.gy=0,0
        global meas_line
        ofile=None
        if format=='txt': #hack for saving in text format
            def save(a,b):
                ouf=open(a+".txt","w")
                if len(b.shape)>1 and min(b.shape)>1:
                    ouf.write("#"+"  ".join(["pnt%i"%i for i in range(1,len(b)+1)])+'\n') #writing header
                savetxt(ouf,b.transpose(),fmt="%8.5f")
                ouf.close()
        elif format=='pts':
            ofile=open(oname,"w")
            ofile.write("[%i,%i]\n"%(self.gx,self.gy))
        #rate=newrate
        if swapaxes: xdir,ydir=2,1
        else: xdir,ydir=1,2
        if control:
            control['size']=(xcnt,ycnt)
            control['nx']=xcnt #is OK for round samples ??
            control['ny']=ycnt
            if 'stop' in control: del control['stop']
        pos_cent=[0,0]
        pos_beg=[self.gx,self.gy]
        if radius>0: # round sample
            if center: self.rate(ydir,-radius)
            if radius<ycnt*ystep-1: ycnt=radius//ystep+1
            if radius<xcnt*xstep-1: xcnt=radius//xstep+1
            xmin,xmax=-xcnt+1,xcnt
            ymin,ymax=-ycnt+1,ycnt
        else:
            xmin,xmax=0,xcnt
            ymin,ymax=0,ycnt
        if testonly: print("scan in range %i-%i/%i-%i"%(xmin,xmax,ymin,ymax))
        nmeas=0
        step_virtual=[0,0]
        for j in range(ymin,ymax):
                meas_line=[]
                pos_line=[]
                k=0
                for i in range(xmin,xmax):
                    if ((radius>0) and ((i*xstep-pos_cent[0])**2+(j*ystep-pos_cent[1])**2)>radius**2) or (centerfirst and i>xmin):
                        meas_line.append(zeros(self.pixtable.shape)) # not measuring at this point - empty data
                        step_virtual[xdir-1]+=xstep
                        continue
                    if k!=0:
                        step_virtual[xdir-1]+=xstep
                        #print('going '+str(xdir)+' '+str(xstep))
                    else: k=1  # first point in the line
                    if step_virtual[xdir-1]!=0: #realize all delayed moves
                        self.rate(xdir,step_virtual[xdir-1])
                    if step_virtual[ydir-1]!=0: #realize all delayed moves
                        self.rate(ydir,step_virtual[ydir-1])
                    step_virtual=[0,0]
                    if control and 'stop' in control: break
                    ############ measurement #######################
                    meas_line.append(self.result())
                    nmeas+=1
                    pos_line.append("[%i,%i]"%(self.gx,self.gy))
                    if recalib>0 and nmeas%recalib==0: # go to calibration point
                        pos_act=[self.gx,self.gy]
                        self.rate(1,pos_calib[0]-pos_act[0])
                        self.rate(2,pos_calib[1]-pos_act[1])
                        #self.rate(xdir,pos_calib[0]-pos_act[0])
                        #self.rate(ydir,pos_calib[1]-pos_act[1])
                        newflat=self.measure()
                        if self.flat!=None:
                            ratio=newflat/self.flat
                            print("recalib. by fact. %.3f"%ratio.mean())
                        self.flat=newflat
                        pos_line.append("[%i,%i]"%(self.gx,self.gy))
                        self.rate(1,-pos_calib[0]+pos_act[0])
                        self.rate(2,-pos_calib[1]+pos_act[1])
                    #print('just measured %i %i'%(i,j))
                    if control:
                        if 'wait' in control: sleep(control['wait'])
                        control['x']=i
                        control['y']=j
                        control['gx']=self.gx
                        control['gy']=self.gy
                        if rc.debug==2: print("at pos %i / %i"%(i,j))
                        elif rc.debug>2: print("at pos %i / %i"%(self.gx,self.gy))
                        if 'meas_evt' in control: control['meas_evt'].set() #starting everything that should come with new measurement
                        if 'queue' in control: control['queue'].put('measured %i %i'%(i,j)) #synchronization for GUI
                        if 'anal_evt' in control: # waiting for analysis to end
                            control['anal_evt'].wait()
                            control['anal_evt'].clear()
                        if 'stop' in control: break
                        if 'refer' in control and control['refer']==1:
                            if j==-ycnt+1: #save calibration
                                if self.flat!=None: self.flat*=array(meas_line[-1])
                                else: self.flat=array(meas_line[-1])
                                self.last=1
                                control['refer']==0
                        if 'anal' in control: control['anal']() #calling analysis directly, not multitasking
                        if 'meas_evt' in control: control['meas_evt'].clear()
                        if rc.single_central_pt and (j==0): break
                if control and 'stop' in control:
                    print('stop forced, leaving')
                    break #leaving
                if j%2==1:
                    meas_line=meas_line[::-1]
                    #print('line %i inverted'%j)
                if ofile!=None:
                    ofile.write("\t".join(pos_line)+"\n")#self.config.Material.decode('ascii')!='simu':
                else:
                    if radius<=0:
                        save(oname+"%03i"%(j+1),array(meas_line))
                    else:
                        save(oname+"%03i"%(j+ycnt),array(meas_line))
                if radius<0: #nevim proc to tu je
                    if j>=radius: break
                    if j<-radius: continue
                    self.rate(xdir,xstep*(int(sqrt(radius**2-(j+1)**2))-int(sqrt(radius**2-j**2))))
                #print("now next line")
                if j<ycnt-1:
                    step_virtual[ydir-1]+=ystep
                if not(rc.polar_stage and rc.single_central_pt and (j==0)):
                    xstep=-xstep
        if control:
            if 'meas_evt' in control:
                control['stop']=2
                control['meas_evt'].set()
            if 'queue' in control: control['queue'].put('finished')
            if 'return' in control: # return to the starting position
                print('returning from pos %i, %i (%i, %i)'%(i,j,xmin,ymin))
                self.rate(ydir,-ystep*(j-ymin))
                ymod=rc.single_central_pt and 1 or 0
                if (j-ymin)%2==ymod: self.rate(xdir,xstep*(i-xmin),wait=-1) #odd number of x-scans
                control['x']=xmin
                control['y']=ymin
                xstep=-xstep
        if hasattr(self,'pixtable') and ofile==None: #saving calibration data
            if self.flat==None: save(oname+"calib",self.pixtable)
            else: save(oname+"calib",array([self.pixtable,self.flat]))
        if returnlast:
            self.rate(1,pos_beg[0]-self.gx)
            self.rate(2,pos_beg[1]-self.gy)
        if ofile!=None:
            if returnlast:  ofile.write("[%i,%i]\n"%(self.gx,self.gy))
            ofile.close()