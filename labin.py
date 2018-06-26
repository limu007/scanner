# -*- coding: utf-8 -*-

testonly=False

######### CRYOSTAT #############
global ser,per
ser=None # serial port connection
per=None # persistent instrument
serspeed=19200 #38400
gloud=0

global meas_line #storage for measurements - made global for testing purposes
meas_line=[]

try:
    import labrc as rc
except:
    from scanner import labrc as rc
    
from numpy import iterable

def unitconv(ufrom,uto,values):
    from scanner.spectra import ev2um
    conv,cinv=1.,0
    if ufrom==uto:
        if not iterable(values): return conv,cinv
        return values
    if ufrom==1: conv,cinv=1/(ev2um*1e-3),1 #nm
    elif ufrom==2: conv,cinv=1/ev2um*1e-4,0 #cm^-1
    if uto==1:
        conv/=ev2um*1e-3
        cinv+=1
    elif uto==2:
        #if cinv%2==1: conv*=1/ev2um*1e-4
        #else:
        conv/=1/ev2um*1e-4
    if not iterable(values): return conv,cinv
    if cinv%2==1: return conv/values
    return values*conv


######### SPECTROSCOPY #############
from ctypes import c_ushort, c_uint8, c_uint16, c_uint32, c_float, c_double, c_char, c_int, pointer
from ctypes import windll, Structure

class DarkCorrectionType(Structure):
    _fields_ = [("m_Enable", c_uint8),
                ("m_ForgetPercentage", c_uint8)]

#class DevStatType(Structure):
#    _fields_ = [("m_SerialId", c_uint16),
#                ("m_UserId", c_uint8),
#                ("m_Stat", c_uint8)]

class IdentType(Structure):
    _fields_ = [("m_SerialId", c_char*10),
                ("m_UserId", c_char*64),
                ("m_Stat", c_uint8)]

class SmoothingType(Structure):
    _fields_ = [("m_SmoothPix", c_uint16),
                ("m_SmoothModel", c_uint8)]

class ControlSettingsType(Structure):
    _fields_ = [("m_StrobeControl", c_uint16),
                ("m_LaserDelay", c_uint32),
                ("m_LaserWidth", c_uint32),
                ("m_LaserWaveLength", c_float),
                ("m_StoreToRam", c_uint16)]

class TriggerType(Structure):
    _fields_ = [("m_Mode", c_uint8),
                ("m_Source", c_uint8),
                ("m_SourceType", c_uint8)]

class MeasConfigType(Structure): #MeasConfigType
    _fields_ = [("m_StartPixel", c_uint16),
                ("m_StopPixel", c_uint16),
                ("m_IntegrationTime", c_float),
                ("m_IntegrationDelay", c_uint32),
                ("m_NrAverages", c_uint32),
                ("m_CorDynDark", DarkCorrectionType),
                ("m_Smoothing", SmoothingType),
                ("m_SaturationDetection", c_uint8),
                ("m_Trigger", TriggerType),
                ("m_Control", ControlSettingsType)]

class SimpleConfigType(Structure): #MeasConfigType
    _fields_ = [("m_StartPixel", c_uint16),
        ("m_StopPixel", c_uint16),
        ("m_IntegrationTime", c_uint16),
        ("m_NrAverages", c_uint32),
        ("m_CorDynDark", DarkCorrectionType),
        ("m_Smoothing", SmoothingType),
        ("Material", c_char*64)]

polycal=[ -3.46472e-10, -5.39531e-06, 3.44907e-01, 1.748616e+02] # fit of pixel2waveleng. dependence for AvaSpec-3648
global specdata
specdata=None
pos_calib=[0,0]

class specscope(object):
    '''generic class what every spectroscope should do
    - by default returns simulated data (for testing purposes)
    '''
    dark=None
    flat=None
    last=None
    step=0.0008267
    pixtable=None
    config=None
    data={}
    pixrange=[0,-1]
    dirs=[0,0]
    #global gx,gy
    gx,gy=0,0
    name=""

    def __init__(self,erange=[1.5,4.5],nbin=None):
        '''range in eV'''
        if nbin==None: nbin=int((erange[1]-erange[0])/self.step)
        else: self.step=float(erange[1]-erange[0])/nbin
        from numpy import arange
        self.pixtable=arange(erange[1],erange[0],-self.step)
        return
    def setup(self,prange=[1,1000],integ=100,aver=10,**kwargs):
        if self.config==None:
            self.config=SimpleConfigType()
        if prange!=None:
            if type(prange[0])==int:
                if prange[1]>=len(self.pixtable):prange[1]=len(self.pixtable)-1
            else:
                prange=[(p>self.pixtable).sum() for p in prange]
                if prange[0]!=0:prange[0]-=1
                #prange[1]=(prange[1]<self.pixtable).sum()
                print('measuring from pix. %i to pix. %i'%tuple(prange))
            if len(self.data)>0: self.data={}
            self.config.m_StartPixel,self.config.m_StopPixel=tuple(prange)
        self.config.m_IntegrationTime=int(integ)
        self.config.m_NrAverages=int(aver)
        self.config.Material=b'simu'
        #if 'aver' in kwargs: self.config.m_NrAverages=int(kwargs['aver'])
        return
    def measure(self,npix=None,**kwargs):
        from numpy.random import normal
        from math import sqrt
        if self.config==None:
            print("You must config first!")
            return
        if 'noise' in kwargs: noise=float(kwargs['noise'])
        else: noise=getattr(self.config,'Noise',0.05)
        noise/=sqrt(self.config.m_NrAverages)
        rep=normal(size=self.pixtable.shape,scale=noise)
        mater=None
        if 'mater' in kwargs: mater=kwargs['mater']
        elif hasattr(self.config,'Material') and self.config.Material.decode('ascii')!='simu': mater=self.config.Material.decode('ascii')
        else: mater=rc.simu_mater
        #print('using |'+str(mater)+'|')
        if mater!=None:
            if type(mater)==list:
                for m in mater:
                    if not m in self.data:
                        if m in rc.eps_ext:
                            from profit import dielect
                            self.data[m]=dielect(self.pixtable,rc.eps_ext[m][0],rc.eps_ext[m][1])
                        elif m in rc.mod_ext:
                            from numpy import polyval
                            self.data[m]=polyval(rc.mod_ext[m],self.pixtable)**2
                        print('calculating diel. '+m)
                if len(mater)-len(rc.simu_layer) in [0,1]:
                    from profit import plate
                    wid=rc.simu_layer
                    if hasattr(rc,'simu_layer_fluc'): wid=[wid[i]*normal(1,rc.simu_layer_fluc[i]) for i in range(len(wid))]
                    rep+=plate(self.pixtable,[self.data[m] for m in mater],wid)
            else:
                if not mater in self.data:
                    if mater in rc.mod_ext:
                        from numpy import polyval
                        self.data[mater]=polyval(rc.mod_ext[mater],self.pixtable)**2
                    else:
                        from scanner.spectra import dbload
                        if type(mater)==bytes: dbname=mater.decode('ascii')
                        else: dbname=mater
                        if dbname in rc.eps_trans: dbname=rc.eps_trans[dbname]
                        try:
                            inp=dbload(dbname,self.pixtable)
                        except:
                            print('cannot load material '+dbname)
                            return rep
                        #for i in range(len(inp)):
                        self.data[mater]=inp[1]
                if len(self.data)>0:
                    from profit import reflect
                    rep+=reflect(self.data[mater])#,self.data[mater[1]])
        rep*=self.config.m_IntegrationTime
        return rep
    def end(self):
        return
    def smooth(self,iata,wid=5):
        from numpy import r_,convolve,hamming
        data=r_[iata[wid//2-1::-1],iata,iata[:-wid//2-1:-1]] #reflection on edges
        box=hamming(wid)
        box/=box.sum()
        return convolve(box,data,"same")[wid//2:-wid//2]

    def result(self,smooth=0,sub_dark=True,div_flat=True,maxit=0):
        '''complete measurement process
        repeats measurement several [maxit] times, until gets positive data
        smooth: applies software smoothing
        dark current and reference correction
        '''
        from numpy import array,r_,iterable
        from time import sleep
        data=None
        iata=None
        if maxit==0: maxit=rc.spectr_nb_trials
        for i in range(maxit):
            iata=self.measure(smooth=smooth)
            if not iterable(iata): continue
            if len(iata)>0:
                data=array(iata)
                if data.max()>0: break
            sleep(0.2)
        if not iterable(data):
                print('measurement failed')
                return # all attempts failed
        #if smooth>0: data=self.smooth(iata,smooth)
        if sub_dark and iterable(self.dark): data-=self.dark
        if div_flat and iterable(self.flat):
                data/=self.flat
                if rc.hard_spec_min: data[data<rc.hard_spec_min]=rc.hard_spec_min
                if rc.hard_spec_max: data[data>rc.hard_spec_max]=rc.hard_spec_max
        self.last=data[:len(self.pixtable)]
        return self.last

    def corr_flat(self,minval=0,minlen=10,alert=None):
        nbin=sum(self.flat<=minval)
        text="%i reference values <b>not positive</b>!"%nbin
        #self.instr.flat=None
        from numpy import int,array,where
        # looking for the longest non-zero interval
        sele=(self.flat>minval).astype(int)
        dele=sele[1:]-sele[:-1]
        p1,p0=list(where(dele>0)[0]),list(where(dele<0)[0])
        p0.insert(0,0)
        p1.append(len(dele))

        if len(p0)<len(p1): p1=p1[:len(p0)]
        elif len(p0)>len(p1): p0=p0[-len(p1):]
        k=(array(p0)-array(p1)).argmax()
        prange=[p1[k]+1,p0[k]+1]
        print("widest nonzero range %i - %i"%tuple(prange))
        if prange[1]-prange[0]>minlen:
            conf=self.instr.config
            conf.m_StartPixel,conf.m_StopPixel=prange[0],prange[1]
            self.instr.pixtable=self.instr.pixtable[prange[0]:prange[1]]
            self.instr.flat=self.instr.flat[prange[0]:prange[1]]
            self.instr.dark=self.instr.dark[prange[0]:prange[1]]
            text+="<br/>reducing spectral range [%.2f - %.2f eV]"%(self.pixtable.min(),self.pixtable.max())
            if alert!=None: alert(self,text)
        else:
            if alert!=None: alert(self,text+"no good calib. range found")
            self.flat[self.flat<=minval]=minval+1.

    def acomm(self,comm,ntries=3,comm_gap=0.1,end_gap=0.05):
        from time import sleep
        for i in range(ntries):
            try:
                self.ard.write(str(comm+"\n"))
            except:
                self.ard.write(bytes(comm+"\n","ascii"))
            sleep(comm_gap)
            try:
                reply=''.join([a.decode("ascii") for a in self.ard.readlines()])
            except:
                reply=''.join(self.ard.readlines())
            if reply.find("nvalid")<0: break
        if gloud>0: print('ARD: '+comm+' [%i]'%i)
        if gloud>1: print('REP: '+reply)
        sleep(end_gap)
        return i

    def adrive(self,peri=10,gap=0):
        import serial
        if not hasattr(self,'ard'):
            try:
                self.ard=serial.Serial(rc.ard_port,baudrate=rc.ard_speed,timeout=2)
            except:
                print("cannot connect to Arduino")
                return -1
        if gap==0: gap=peri//2
        if rc.ard_mode==1:
            rep=self.acomm('SING')
            if rep>2: print("Arduino init failed")
            else: print("single drive mode")
        if peri: self.acomm('PERI %i'%peri)
        if gap: self.acomm('GAP %i'%gap)
        rc.motor_lin_speed=rc.motor_lin_speed*peri/5.
        print("high %i ms low %i ms"%(peri-gap,gap))
        self.ard.readlines()
        return 0

    def rate(self,ch=1,nstep=100,wait=0,ain=0,ntries=3,loud=0):
        '''controlling scanner
        ch: defines movement axis
        adapted for different drives:
            Thorlabs rotation stage
            A-drive linear motor controlled by Avantes
                                            by Arduino
        '''
        from time import sleep
        #global gx,gy
        if hasattr(self,'motor') and ch==2: # rotation stage by Thorlabs
            self.gy+=nstep
            if testonly: return
            self.motor.relat(nstep)
            if wait>=0: sleep(rc.motor_rot_wait+rc.motor_rot_speed*abs(nstep)) # needs some time to reach the goal
            # for the moment we are not able to check the time of arrival
            # should probably use GetPosition function
            return 0
        if hasattr(self,'ard') and self.ard!=None: # control by Arduino
            if wait<=0: wait=rc.comm_ard_wait
            if ch==1: self.gx+=nstep
            else: self.gy+=nstep
            if testonly: return
            dirtext=""
            if nstep<0:
                nstep=-nstep
                if self.dirs[ch-1]!=-1:
                    if ch==2: dirtext="LEFT"
                    else: dirtext="UP"
                    self.dirs[ch-1]=-1
            else:
                if self.dirs[ch-1]!=1:
                    if ch==2: dirtext="RIGHT"
                    else: dirtext="DOWN"
                    self.dirs[ch-1]=1
            if len(dirtext)>0:
                self.acomm(dirtext)
                sleep(wait)
            step_so_far=0
            for i in range(ntries): # X atempts to communicate
                text="STEP%i %i \r\n"%(ch,nstep-step_so_far)
                if rc.ard_comver==3: text=bytes(text,"ascii")
                self.ard.write(text)
                if gloud>0: print("STEP%i %i [ %i %s]"%(ch,nstep-step_so_far,self.dirs[ch-1],dirtext))
                wait=rc.motor_lin_wait+rc.motor_lin_speed*abs(nstep)
                sleep(wait)
                #sleep(wait)
                if rc.ard_comver==3: reply=''.join([a.decode("ascii") for a in self.ard.readlines()])
                else: reply=''.join(self.ard.readlines())
                if not reply.find("finished")>0:
                    for j in range(5):
                        sleep(rc.motor_lin_wait)
                        if rc.ard_comver==3: reply+=''.join([a.decode("ascii") for a in self.ard.readlines()])
                        else: reply+=''.join(self.ard.readlines())
                        if reply.find("finished")>0: break
                if reply.find("nvalid")>0 or loud>2: print(" || ".join(reply.split("\n")))
                if reply.find("nvalid")>0:
                    p0=reply.find("argv[0] =>")
                    if p0>0 or gloud>1:
                        print("analyzing:",reply[p0+7:p0+17])
                    p1=reply.find("argv[1] =>")

                #if reply.find("finished")>0: break
                else: break
                if loud>1:
                    print("retrying\n")
                    print("["+reply+"]")
            return i
        return -1

    def write(self,fname,expername=None, comm=[],data=None,xdata=None,format='%f',type='ascii'):
        '''saving data in ascii format'''
        if data==None: data=self.last
        colab='#'
        if type=='matlab': colab='%%'
        file=open(fname,"w")
        if expername: file.write(colab+" %s\n"%expername)
        #if expernote!='': file.write(comm+" %s\n"%expernote)
        #file.write(comm+" INTEG %.1f AVG %.1f SMOOTH %i\n"%(self.config.m_IntegrationTime,self.mean,self.smooth))
        for li in comm: file.write(colab+li)
        from numpy import iterable
        if xdata==None:xdata=self.pixtable
        multicol=iterable(data[0]) #are data 1-D or more
        for i in range(len(data),0,-1):
            if type=='matlab':
                if multicol>0: lin=(format+rc.output_separ+"%f"+rc.output_separ+"%s\n")%(xdata[i-1],data[i-1][-1],rc.output_separ.join([format%d for d in data[i-1][:-1]]))
                else: lin=(format+rc.output_separ+"%f"+"\n")%(xdata[i-1],data[i-1])
            else:
                if multicol>0: lin=("%f"+rc.output_separ+"%s\n")%(xdata[i-1],rc.output_separ.join([format%d for d in data[i-1]]))
                else: lin=("%f"+rc.output_separ+format+"\n")%(xdata[i-1],data[i-1])
            file.write(lin)
        file.close()

    def pulse(self,peri=10,dur=100,oname='C:\\Users\\Lenovo\\Documents\\Data\\temp_line',loud=0):
        from datetime import datetime
        from time import sleep
        start=datetime.now()
        dela=peri/2.
        act=datetime.now()
        dis=(act-start).total_seconds()
        cnt=0
        while dis<dur:
            if loud:
                sys.stdout.write("Tick %f: %f\n"%(dis,dela))
                sys.stdout.flush()
            sleep(dela)
            act=datetime.now()
            dis=(act-start).total_seconds()
            if dis//peri>cnt:
                vals=self.result()
                cnt+=1
                dela=peri/2.
                if oname: self.write(oname+"_%03i.dat"%cnt,comm=[str(act)+"\r\n"])#,array([self.pixtable,self.last]))
                print("time count %i"%cnt)
                if self.parent: self.parent.graph.update_measure(instr=self,vals=vals)
            else:
                dela=(dis//peri+1)*peri-dis-0.5
                if dela<0: dela=0.1

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

class ocean(specscope):
    '''developped and tested for JAZ/NIRquest
    '''
    hnd=None
    device=None
    #dark=None
    #flat=None
    #last=None
    time=None
    #ddev=None
    #alldev=None
    dsize=2048
    def __init__(self):
        import os,jpype
        if rc.java_jdk: #linux?
            os.environ['OOI_HOME']=rc.ooihome
            os.environ['JAVA_HOME']=rc.java_home
            os.environ['JVM_ROOT']=os.environ['JAVA_HOME']
            os.environ['JDK_INCLUDE_FILE_ROOT']=rc.java_jdk
            os.environ['LD_LIBRARY_PATH']=os.environ['OOI_HOME']+':'+os.environ['JAVA_HOME']+'/amd64/server'

        jpype.startJVM(rc.java_jvm,"-Djava.class.path="+os.environ['OOI_HOME']+"/OmniDriver.jar")
        self.device = jpype.JPackage("com").oceanoptics.omnidriver.api.wrapper.Wrapper();
        self.device.openAllSpectrometers()
        self.dsize=self.device.getNumberOfPixels(0)
        self.name=self.device.getFirmwareModel(0)

    def setup(self,prange=[1,3000],integ=50,aver=10,calib=None,smooth=None,dyndark=False,unit=0):
        self.device.setScansToAverage(0, aver)
        self.device.setIntegrationTime(0, 1000*integ);
        from numpy import polyval,r_

        #rep=self.device.AVS_GetNumPixels(self.hnd,pointer(inum))
        if self.pixtable==None:
            coef=list(self.device.getCalibrationCoefficientsFromEEProm(0).getWlCoefficients())
            sdev=polyval(coef[::-1],r_[:self.dsize])
            self.pixtable=unitconv(1,unit,r_[list(sdev)])
        if dyndark: self.device.setCorrectForElectricalDark(0,rc.dyndark_forget)
        if rc.cool_temp<0:
            twe=self.device.getFeatureControllerThermoElectric(0)
            twe.setDetectorSetPointCelsius(rc.cool_temp)

    def measure(self):
        from numpy import array
        spect = array(self.device.getSpectrum(0))
        return list(spect[:self.dsize])

    def calibrate(nchan=0,poly=[177.7322,0.380633,-0.00001346845,2.9798580E-9]):
        c=self.device.CreateCoefficients()
        c.setWlIntercept(poly[0])
        c.setWlFirst(poly[1])
        c.setWlSecond(poly[2])
        c.setWlThird(poly[3])
        c.setStrayLight(3.0)
        self.device.insertKey("Mat429sky")
        self.device.setCalibrationCoefficientsIntoEEProm(0,nchan,c,True, True, True)

    def end(self):
        import os,jpype,time
        self.device.closeAllSpectrometers()
        #time.sleep(2)
        #del self.device
        jpype.shutdownJVM()

class oceanjaz(ocean):
    '''multichannels
    '''
    nchan=1
    chord=[]
    chrange=[]

    def setup(self,prange=[1,3000],integ=50,aver=10,calib=None,smooth=None,dyndark=False,unit=0):
        from numpy import polyval,r_,array,argsort
        ext=self.device.getWrapperExtensions()
        self.nchan=ext.getNumberOfEnabledChannels(0)
        gchan=[]
        coefbeg=[self.device.getWavelengthIntercept(0,i) for i in range(self.nchan)]
        self.chrange=[[None,None] for i in range(self.nchan)]
        self.chord=list(argsort(coefbeg).astype(int))
        print
        iprev=None
        for i in self.chord:
            coef=array(self.device.getCalibrationCoefficientsFromEEProm(0,int(i)).getWlCoefficients())
            xchan=list(polyval(coef[::-1],r_[:self.dsize]))
            print("ch.%i:%.1f-%.1f"%(i,xchan[0],xchan[-1]))
            if len(gchan)>0:
                nleft=sum(xchan<gchan[-1])
                xmid=nleft//2
                nright=sum(xchan[xmid]<gchan)
                self.chrange[iprev][1]-=nright
                gchan=gchan[:-nright]+xchan[xmid:]
            else:
                xmid=0
                gchan=xchan

            self.chrange[i]=[xmid,self.dsize]
            #print("chan. %i:%.2f-%.2f nm"%(i+1,gchan[-1][0],gchan[-1][-1]))
            self.device.setIntegrationTime(0, int(i), 1000*integ)
            self.device.setScansToAverage(0, int(i), aver)
            if dyndark: self.device.setCorrectForElectricalDark(0,int(i),rc.dyndark_forget)
            iprev=i
        self.pixtable=unitconv(1,unit,array(list(gchan)))

    def measure(self):
        data=[]
        from numpy import array
        for i in self.chord:
            #spectrum = array(self.device.getSpectrum(0,i))
            data+=list(self.device.getSpectrum(0,int(i)))[self.chrange[i][0]:self.chrange[i][1]]
            #print("chan. %i:max %.2f  mean %.2f"%(i+1,spectrum.max(),spectrum.mean()))
        return data

class avantes(specscope):
    '''developped and tested for AvaSpec 3648
    '''
    hnd=None
    device=None
    #dark=None
    #flat=None
    #last=None
    time=None
    ddev=None
    alldev=None
    def __init__(self,port=0,loud=0):
        '''port=0: USB, others are COM
        '''
        self.device = windll.as5216
        l_Port = self.device.AVS_Init(port)
        if loud>0: print("testing device at port %i"%l_Port)
        ndev=self.device.AVS_GetNrOfDevices()
        if ndev<1: return
        fdev = c_double(ndev)
        idev = (c_int*ndev)()
        types=(IdentType*ndev)() #types= (c_int*ndev*100)()
        rep=self.device.AVS_GetList(2*75*ndev,idev,types)
        self.alldev=types
        if rep==0: print('cannot get self.device list')
        elif loud==1: print('found device '+types[0].m_UserId)
        self.hnd=self.device.AVS_Activate(pointer(types[0]))

    def setup(self,prange=[1,3000],integ=50,aver=10,calib=None,smooth=None,dyndark=False,unit=0):
        inum=c_ushort()
        starting=True
        rep=self.device.AVS_GetNumPixels(self.hnd,pointer(inum))
        if rep<0 or inum.value<=10:
            print('invalid number of pixels')
            return rep
        from scanner.spectra import ev2um
        if calib==None: calib=rc.spec_pixtable_source
        if self.config==None:
            starting=True
            self.config=MeasConfigType()
        from numpy import arange,array,polyval
        if prange!=None:
            if calib==1:
                sdev = (c_double*inum.value)()
                self.device.AVS_GetLambda(self.hnd,pointer(sdev))
                self.pixtable=unitconv(1,unit,array(list(sdev)))
            else:
                self.pixtable=unitconv(1,unit,polyval(polycal,arange(inum.value)))
            print("converting pixels into %i"%(unit))
            if type(prange[0])!=int: # float - in physical units
                prange=[(p<self.pixtable).sum() for p in prange[::-1]]
                if prange[0]!=0:prange[0]-=1
                #prange[1]=(prange[1]<self.pixtable).sum()
                print('measuring from %i to %i'%tuple(prange))
            self.pixrange=prange
            self.ddev=None
            self.config.m_StartPixel=prange[0]
            self.config.m_StopPixel=prange[1]
            self.pixtable=self.pixtable[prange[0]:prange[1]]
            self.dark=None
            self.flat=None
        self.config.m_NrAverages=aver
        if False:  #short test measurement
            self.config.m_IntegrationTime=1
            mok=self.device.AVS_PrepareMeasure(self.hnd,pointer(self.config))
            if mok<0:
                print('first setup failed')
            else:
                self.measure()
                print('test measurement passed')
        self.config.m_IntegrationTime=integ
        if dyndark:
            self.config.m_CorDynDark.m_Enable=1
            self.config.m_CorDynDark.m_ForgetPercentage=rc.dyndark_forget
        else:
            self.config.m_CorDynDark.m_Enable=0
        if smooth!=None: self.config.m_Smoothing.m_SmoothPix=smooth
        mok=self.device.AVS_PrepareMeasure(self.hnd,pointer(self.config))
        if mok<0:
            print('setup failed')
            return mok
        self.measure()
        return inum.value

    def measure(self,npix=None,winhand=0,nmeas=1,prep_task=None):
        from time import sleep
        if prep_task!=None: eval(prep_task)
        if winhand<0 and hasattr(self,"parent"): winhand=self.parent
        zoo=self.device.AVS_Measure(self.hnd,winhand,nmeas)
        if zoo==-21:
            print('invalid state for spectra acquisition')
            return
        sleep(rc.sleepcon+rc.sleepfact*self.config.m_IntegrationTime*self.config.m_NrAverages)
        if npix==None: npix=self.config.m_StopPixel-self.config.m_StartPixel
        if self.time == None: self.time = c_int()
        if self.ddev == None: self.ddev = (c_double*(npix+1))()
        out=self.device.AVS_GetScopeData(self.hnd,pointer(self.time),pointer(self.ddev))
        if out==-8: print('no data measured')
        elif out==-4: print('bad handle')
        elif out<0: print('measurement failed [%i]'%int(out))
        else: return list(self.ddev[:npix])

    def revert(self,data=None,newbins=None):
        # gets data for uniform sampling
        from scipy import interpolate
        from numpy import arange
        if data==None: data=self.last
        if self.pixtable[0]<self.pixtable[-1]:
            bins=self.pixtable
        else:
            bins=self.pixtable[::-1]
            data=data[::-1]
        tck=interpolate.splrep(bins,data)
        if newbins==None:
            step=(bins[-1]-bins[0])/(len(bins)-1)
            newbins=arange(bins[0],bins[-1]+step,step)
        return interpolate.splev(newbins,tck)

    def end(self):
        self.device.AVS_StopMeasure(self.hnd)
        self.device.AVS_Deactivate(self.hnd)
        self.device.AVS_Done()

    def do_all(self):
        if self.hnd==None: self.__init__()
        rep=setup(aver=50)
        if rep!=None:
            from numpy import array
            return array(self.measure())
    def rate(self,ch=1,nstep=100,wait=0,ain=0,ntries=3,loud=0):
        if super(avantes, self).rate(ch,nstep,wait,ain,ntries,loud)>=0: return
        # control by Avantes ports
        dir=0
        if type(nstep)==float: # if not converted before
                nstep=int(nstep/rc.xstepsize)
                #print('going %i steps'%nstep)
        if ch==1: self.gx+=nstep
        else:
            self.gy+=nstep
            ch=2
        if nstep<0:
            self.device.AVS_SetDigOut(self.hnd,ch,1)
            dir=1
            nstep=-nstep
        else: self.device.AVS_SetDigOut(self.hnd,ch,0)
        meas=[]
        if ain>0: dd=c_float()

        for i in range(nstep):
                self.device.AVS_SetDigOut(self.hnd,ch+3,1)
                if wait>0: sleep(wait)
                if ain>0:
                        self.device.AVS_GetAnalogIn(self.hnd,ain,pointer(dd))
                        meas.append(dd)
                self.device.AVS_SetDigOut(self.hnd,ch+3,0)
                if wait>0: sleep(wait)
        if dir==1:
                self.device.AVS_SetDigOut(self.hnd,ch,0)
        if ain>0:
                from numpy import array
                return array(meas)

#-----------------------------------------------------------

class linuscope(avantes):
    tmpfile="/tmp/data.txt"
    logfile="/tmp/usblog.txt"
    loglen=0
    loghand=None
    pixorig=None
    repeat=1

    def __init__(self,comm="usbdir",loud=0):
        '''communicating through C++ interface on linux
        '''
        import os
        import subprocess as sp
        if hasattr(rc,'scope_pixtable'): pixfile=rc.scope_pixtable
            #default pixel position table (under linux, we don't receive pixel positions)
        else:
            pixfile=self.tmpfile
        import time
        is_ok=False
        if self.logfile:
            if os.path.exists(self.logfile): os.unlink(self.logfile)
            self.loghand=open(self.logfile,"w")
        else:
            self.loghand=sp.PIPE
        self.process = sp.Popen(comm,stdin=sp.PIPE,stdout=self.loghand,universal_newlines=True)#os.popen(comm,"w")
        self.device = self.process.stdin
        for j in range(5):
            self.device.write("MEAS %s\n"%self.tmpfile);self.device.flush()
            for i in range(10):
                time.sleep(0.5)
                if os.path.exists(self.tmpfile):
                    is_ok=True
                    break
            if is_ok and os.path.getsize(self.tmpfile)>1000: break
            else: is_ok=False
        if not is_ok:
            print('init failed')
            self.pixtable=None
            return
        data=[a.strip().split() for a in open(pixfile).readlines() if a[0]!='#']
        if len(data)<=1: #no measurement succeeded
            self.pixtable=None
            return
        from numpy import array
        from scanner.spectra import ev2um
        self.pixorig=1000./ev2um/array([float(a[0]) for a in data])
        self.pixtable=self.pixorig
        self.loghand.flush()
        self.loglen=len(open(self.logfile).readlines())

    def setup(self,prange=[1,1000],integ=100,aver=10,calib=0,smooth=None,dyndark=None,unit=0):
        if self.config==None: self.config=MeasConfigType()
        if prange!=None:
            if type(prange[0])==int: # really in pixels
                if prange[1]>=len(self.pixorig):prange[1]=len(self.pixorig)-1
            else: # in eV?
                #prange=[(p>self.pixorig).sum() for p in prange]
                from numpy import where
                print("(got to convert "+str(prange)+")")
                irange=[where(p<self.pixorig)[0] for p in prange]
                prange=[(len(p)>0 and p[-1] or 0) for p in irange]
            if prange[0]>prange[1]:prange=prange[::-1]
            if prange[0]>0: prange[0]-=1
            if abs(prange[1]-prange[0])==0: prange=[0,len(self.pixorig)-1]
            #prange[1]=(prange[1]<self.pixtable).sum()
            print('measuring from pix. %i to pix. %i'%tuple(prange))
            if len(self.data)>0: self.data={}
            self.pixrange=prange
            self.config.m_StartPixel=prange[0]
            self.config.m_StopPixel=prange[1]
            self.pixtable=self.pixorig[prange[0]:prange[1]] #convert to right units?
            self.pixtable=unitconv(0,unit,self.pixtable)
            self.dark=None
            self.flat=None
        if dyndark!=None:
            if dyndark==True: dyndark=50
            if type(dyndark)==int and dyndark>0:
                self.config.m_CorDynDark.m_Enable=1
                self.config.m_CorDynDark.m_ForgetPercentage=dyndark
                self.device.write("DYD %i\n"%50);self.device.flush();
            else:
                self.config.m_CorDynDark.m_Enable=0
                self.device.write("DYD %i\n"%0);self.device.flush();
        if smooth!=None:
            self.config.m_Smoothing.m_SmoothPix=smooth
            self.device.write("SMO %i\n"%smooth);self.device.flush();
        self.config.m_IntegrationTime=int(integ)
        if aver>0:
            if rc.max_aver_count and aver>rc.max_aver_count:
                self.repeat=2*aver//rc.max_aver_count
                aver=rc.max_aver_count//2
                print("repeating %ix"%self.repeat)
            else:
                self.repeat=1
            self.config.m_NrAverages=aver
            self.device.write("AVG %i\n"%self.config.m_NrAverages);self.device.flush();
        if integ>0:
            self.config.m_IntegrationTime=integ
            self.device.write("EXP %i\n"%integ);self.device.flush();

        return


    def measure(self,npix=None,**kwargs):
        import os,time
        if os.path.exists(self.tmpfile): otime=os.path.getmtime(self.tmpfile)
        else: otime=0
        from numpy import array,zeros
        if self.pixrange[1]<0: dsum=zeros(len(self.pixtable))
        else: dsum=zeros(self.pixrange[1]-self.pixrange[0])
        for i in range(self.repeat):
            self.device.write("MEAS %s\n"%self.tmpfile);self.device.flush()
            if self.repeat>1:
                if self.parent and hasattr(self.parent,"pbar"): self.parent.pbar.setValue(((i+1)*100)/self.repeat)
                #print("it %i:"%(i+1),)
            if self.config: time.sleep((self.config.m_IntegrationTime-1)/1000.)
            for i in range(30):
                self.loghand.flush()
                actlog=open(self.logfile).readlines()
                if self.loglen<len(actlog):
                    self.loglen=len(actlog)
                    if actlog[-1].find('closing')>0 or actlog[-2].find('closing')>0: break
                #if os.path.exists(self.tmpfile) and otime<os.path.getmtime(self.tmpfile): break
                time.sleep(0.3)
                #!!!! need to check the creation time
            time.sleep(rc.linux_write_delay)
            if not os.path.exists(self.tmpfile):
                print("file %s missing"%self.tmpfile)
                return None
            for i in range(3):
                time.sleep(rc.linux_write_delay)
                data=[a.strip().split() for a in open(self.tmpfile).readlines() if a[0]!='#']
                if len(data)==len(self.pixorig): break
            dsum+=array([float(a[-1]) for a in data if len(a)==2][self.pixrange[0]:self.pixrange[1]])
            self.ddev=dsum/self.repeat
        return self.ddev

    def end(self):
        self.device.write("END\n");#self.device.flush()
        self.device.close()

    def check_bord(self,sat_lev=14000):
        #border reached
        from numpy import array
        data=self.last
        if sum(data>sat_lev)>10:
            print('saturating')
            self.config.m_IntegrationTime=int(self.config.m_IntegrationTime*0.8)

def webfetch(url,timeout=4,ver=1):
    if ver==1:
        from urllib import request
        try:
            http=request.urlopen(url, timeout=timeout)
        except:
            return None
        return http.read().decode('utf-8')
    import urllib3
    http = urllib3.PoolManager()
    try:
        link=http.request('GET',url,pool_timeout=timeout,retries=1)
    except:
        return None
    return link.data

def overlaps(caltab):
    from numpy import sum,where
    cn=(sum([c.sum(0) for c in caltab],0)-0.9).astype(int)
    istart=where(cn[1:]>cn[:-1])[0]
    iend=where(cn[1:]<cn[:-1])[0]
    return list(istart),list(iend)

def gettrans(enx,dpos,ref,smot=0.03,skiplowest=0,rep=1,osel=None,disfun=None,weifun=None):
    ysel=None
    from numpy import array,percentile,iterable,newaxis
    if skiplowest>0:
        ysel=[ref[i]>percentile(ref[i],skiplowest) for i in range(len(ref))]
        if iterable(osel) and len(osel)>=len(ref):
            for i in range(len(ref)):
                ysel[i]*=osel[i]
        ref=[ref[i][ysel[i]] for i in range(len(ref))]
    else:
        if iterable(osel) and len(osel)>=len(ref):
            ysel=osel
    caltab=[]
    for j in range(len(ref)):
        if iterable(ysel):
            cmat=smot-abs(dpos[j][ysel[j]][:,newaxis]-enx[newaxis,:])
        else:
            cmat=smot-abs(dpos[j][:,newaxis]-enx[newaxis,:])
        cmat[cmat<0]=0
        if disfun!=None:
            cmat=disfun(cmat)
        csum=cmat.sum(0)
        cmat[:,csum>0]/=csum[csum>0]
        caltab.append(cmat)
    if rep==-1: return caltab
    istart,iend=overlaps(caltab)
    for i in range(len(iend)):
        print(iend[i]-istart[i])
    for i in range(len(istart)):
        zub=array([ref[j].dot(caltab[j])[istart[i]:iend[i]+1] for j in range(len(ref))])
        if weifun!=None:
            for k in range(len(zub)):
                if sum(zub[k])==0: continue
                zub[k]=weifun(zub[k])
        print(zub.shape)
        zn=zub.sum(1)
        zn[zn==0]=1
        zub/=zn[:,newaxis]
        zn=zub.sum(0)
        zub/=zn[newaxis,:]
        for j in range(len(caltab)):
            caltab[j][:,istart[i]:iend[i]+1]*=zub[j]
    return caltab,ysel

def gettransold(enx,dpos,ref,smot=0.03,skiplowest=0,rep=1,scale=[],spow=1):
    '''create transformation table
    '''
    import numpy as np
    ysel=None
    if skiplowest>0:
        ysel=[ref[i]>np.percentile(ref[i],skiplowest) for i in range(len(ref))]
        ref=[ref[i][ysel[i]] for i in range(len(ref))]
    else:
        ysel=[d>0 for d in dpos]
    caltab=[]
    for j in range(len(ref)):
        cmat=smot-abs(dpos[j][ysel[j]][:,np.newaxis]-enx[np.newaxis,:])
        cmat[cmat<0]=0
        csum=cmat.sum(0)
        cmat[:,csum>0]/=csum[csum>0]
        caltab.append(cmat)
    if rep==-1: return caltab,ysel
    istart,iend=overlaps(caltab)
    for i in range(len(iend)):
        print(iend[i]-istart[i])
    for k in range(len(istart)):
        zub=np.array([ref[j].dot(caltab[j])[istart[k]:iend[k]+1] for j in range(len(ref))])
        #if len(scale)==len(zub): zub=np.array(scale)[:,np.newaxis]*zub
        for i in range(len(zub)):
            zub[i]-=zub[i].min()
            if zub[i].max()==0: continue
            if spow<0:
                minpos=np.argmin(zub[i])
                maxpos=np.argmax(zub[i])
                print("[%i:%i]"%(minpos,maxpos))
                if minpos<maxpos:
                    if minpos>0: zub[i][:minpos]=0
                    if maxpos<len(zub[i])-1: zub[i][maxpos+1:]=zub[i][maxpos]
                else:
                    if minpos<len(zub[i])-1: zub[i][minpos:]=0
                    if maxpos>0: zub[i][:maxpos]=zub[i][maxpos]
        #print(zub.shape)
        zn=zub.sum(1)
        zn[zn==0]=1
        zub/=zn[:,np.newaxis]
        if spow>1:
            zub=zub**abs(spow)
        zn=zub.sum(0)
        zn[zn==0]=1
        zub/=zn[np.newaxis,:]
        for j in range(len(caltab)):
            caltab[j][:,istart[k]:iend[k]+1]*=zub[j]
    if len(scale)==len(caltab):
        for j in range(len(caltab)):
            caltab[j]/=scale[j]
    return caltab,ysel

loud=0
class webocean(specscope):
    path="http://localhost:8889/?exp="
    schan=0
    intime=100
    scomm=None
    normfact=[]
    normint=[]
    chanrange=[]
    chanpos=[]
    chanene=[]
    intfact=[]

    darkcorr=False
    nonlincorr=False
    
    def start(self,path,incr=0):
        self.path=path
        comm='java -cp ".;'+rc.java_ooijar+'"  '
        portnum=self.path[self.path.rfind(":")+1:]
        portnum=portnum[:portnum.find('/')]
        if incr!=0:
            portnew=str(int(portnum)+incr)
            self.path=self.path.replace(portnum,portnew)
            portnum=portnew
        comm+='specserve2 '+portnum#.specserve2'
        print("starting "+comm)
        import subprocess,os
        os.environ['PATH']+=";"+rc.jdkpath+";"+rc.java_ooipath
        self.scomm=subprocess.Popen(comm,cwd=rc.java_runpath,stdout=subprocess.PIPE)
        from scanner.exutils import NonBlockingStreamReader as nbsr
        self.spout=nbsr(self.scomm.stdout)
        self.splog=[]
        

    def __init__(self,path=None,loud=0,chan=None):
        self.config=MeasConfigType()
        if path!=None: self.path=path
        if chan!=None: self.schan=chan
        #out=self.acomm.stdout.read()
        self.start(self.path)
        from time import sleep
        sleep(5)
        #print(out)

    def getresp(self,path,minlen=10):
        dstring=webfetch(path,timeout=rc.web_timeout)
        import xml.etree.ElementTree as etree
        if dstring==None:
            print("trying to restart webservice at "+path)
            self.scomm.terminate()
            from time import sleep
            sleep(rc.web_restart_gap)
            oldpath=self.path
            self.start(self.path,incr=1)
            path=path.replace(oldpath,self.path)
            #print("query from "+path)
            dstring=webfetch(path,timeout=rc.web_timeout)
            if dstring==None:
                print("failed")
                return
        #dstring=data.read()
        while True:
            gline=self.spout.readline(0.1)
            if gline==None: break
            self.splog.append(gline.decode())
        if len(self.splog)>2000: self.splog=self.splog[-1000:]
        if len(dstring)<minlen:
            print("truncated data: exiting\n"+dstring.decode('u8'))
            return
        if loud>1: print(dstring.decode('u8'))
        try:
            rep=etree.fromstring(dstring)
            return rep
        except:
            print("XML not valid")
            return dstring

    def setup(self,prange=[1,1000],integ=100,aver=10,calib=0,smooth=None,dyndark=None,unit=0):
        #import urllib2
        from numpy import array,zeros,ones
        tree=self.getresp(self.path+"0&chan=%i"%(self.schan))
        self.intime=integ
        pixval={}
        if tree==None:
            print("spectrometer not working")
            return
        if self.schan>=0:
            self.pixtable=array([float(a) for a in tree.find("data").text.replace(",",".").split()])
            self.pixtable=unitconv(1,unit,self.pixtable)
        else: #all channels
            for ch in tree.findall("data"):
                if not "chan" in ch.attrib: continue
                pixval[int(ch.attrib["chan"])]=array([float(a) for a in ch.text.replace(",",".").split()])
            chanext=zeros((len(pixval),2))
            for i in range(len(pixval)):
                chanext[i]=pixval[i][[0,-1]]
            if chanext[0][0]>chanext[0][1]: chanext=chanext[:,::-1]
            sord=chanext[:,0].argsort()
            print(chanext[sord])
            chsep=(chanext[sord][:-1,1]+chanext[sord][1:,0])/2.
            urange=unitconv(unit,1,array(prange))
            urange.sort()
            print(urange)
            chsep=[urange[0]]+list(chsep)+[urange[1]]
            tlen=0
            print(chsep)
            self.chanrange=[[]]*len(sord)
            for i in range(len(sord)):
                imin=sum(chsep[i]>=pixval[sord[i]])
                imax=sum(chsep[i+1]>=pixval[sord[i]])
                tlen+=imax-imin
                self.chanrange[sord[i]]=[imin,imax]
            print([pixval[i][[self.chanrange[i][0],self.chanrange[i][1]-1]] for i in range(len(sord))])
            self.pixtable=zeros(tlen)
            pint=0
            self.chanpos=[[]]*len(sord)
            self.chanene=[unitconv(1,unit,pixval[i]) for i in range(len(sord))]
            self.chanval=[[]]*len(sord)
            for i in sord:
                dint=self.chanrange[i]
                self.chanpos[i]=[pint,pint+dint[1]-dint[0]]
                pint+=dint[1]-dint[0]
            for i in sord:
                dint=self.chanrange[i]
                pint=self.chanpos[i]
                self.pixtable[pint[0]:pint[1]]=pixval[i][dint[0]:dint[1]]
            self.normfact=ones(len(sord))
            self.pixtable=unitconv(1,unit,self.pixtable)
        if type(prange[0])!=int: # float - in physical units
            print("converting range %.2f - %.2f"%tuple(prange))
            prange=[(p<self.pixtable).sum() for p in prange[::-1]]
            if prange[0]!=0:prange[0]-=1
            if prange[1]<prange[0]: prange=prange[::-1]
            self.config.m_StartPixel,self.config.m_StopPixel=prange
            print("chosen range %i:%i"%tuple(prange))
        self.config.m_NrAverages=int(aver)
        return pixval
        
    def makequery(self,chan=-1):
        if chan<0: chan=self.schan
        if len(self.intfact)>0:
            query=self.path.replace("exp=","&".join(["exp%i=%i"%(i,int(self.intfact[i]*self.intime)) for i in range(len(self.intfact))]))
            query+="&chan=%i"%chan
        else:
            query=self.path+"%i&chan=%i"%(self.intime,chan)
        query+="&corr=%i"%(2*self.darkcorr+self.nonlincorr)

        if self.config.m_NrAverages>1: query+="&avg=%i"%self.config.m_NrAverages
        if loud>0: print(query)
        return query
        
    def measure(self,npix=None,**kwargs):
        #import urllib2
        import xml.etree.ElementTree as etree
        from numpy import array,zeros_like,iterable
        
        #data=urllib2.urlopen(query)
        tree=self.getresp(self.makequery())
        if tree==None: return 
        if self.schan>=0:
            self.ddev=array([float(a) for a in tree.find("data").text.replace(",",".").split()])
            return self.ddev[self.config.m_StartPixel:self.config.m_StopPixel]
        self.ddev=zeros_like(self.pixtable)
        for ch in tree.findall("data"):
            if not "chan" in ch.attrib: continue
            id=int(ch.attrib["chan"])
            if id>=len(self.chanrange): continue
            dint=self.chanrange[id]
            pint=self.chanpos[id]
            if pint[0]>=pint[1]: continue #disabled channel
            pixval=array([float(a) for a in ch.text.replace(",",".").split()])
            self.chanval[id]=pixval
            if "smooth" in kwargs and kwargs["smooth"]>0: pixval=self.smooth(pixval,kwargs["smooth"])
            self.ddev[pint[0]:pint[1]]=pixval[dint[0]:dint[1]]*self.normfact[id]
        return self.ddev

    def saveall(self,basename):
        from numpy import savetxt,array
        for i in range(len(self.chanpos)):
            savetxt(basename+"_%02i.txt"%i,array([self.chanene[i],self.chanval[i]]).T, fmt="%8.5f")
            
    def end(self):
        if self.scomm!=None:
            del self.spout 
            self.scomm.terminate()
            #self.scomm.send_signal(15)

class uniocean(webocean):
    transtable=None
    ysel=None
    
    def __init__(self,path=None,loud=0,chan=None):
        self.darkcorr=True
        self.nonlincorr=True
        webocean.__init__(self,path,loud,chan)
        
    def setup(self,prange,estep=0.02,integ=100,aver=10,unit=0,minpix=10):
        from numpy import array,zeros,ones,arange
        self.pixtable=arange(prange[0],prange[1],estep)
        tree=self.getresp(self.path+"0&chan=%i"%(self.schan))
        self.intime=integ
        pixval={}
        if tree==None:
            print("invalid calib. response")
            return -1
        for ch in tree.findall("data"):
            if not "chan" in ch.attrib: continue
            pixval[int(ch.attrib["chan"])]=array([float(a) for a in ch.text.replace(",",".").split()])
        self.intfact=list(ones((len(pixval))))
        chanext=zeros((len(pixval),2))
        for i in range(len(pixval)):
            chanext[i]=pixval[i][[0,-1]]
            #if sum(chanext[i])>prange[1]
        if chanext[0][0]>chanext[0][1]: chanext=chanext[:,::-1]
        sord=chanext[:,0].argsort()
        self.chanene=[unitconv(1,unit,pixval[i]) for i in range(len(sord))]
        self.chanval=zeros((len(sord),max([len(c) for c in self.chanene])))
        for i in range(len(pixval)):
            if sum(self.chanene[i]>prange[0])<minpix: self.intfact[i]=0
            if sum(self.chanene[i]<prange[1])<minpix: self.intfact[i]=0
        self.config.m_NrAverages=int(aver)

    def measure(self):
        #self.ysel=None
        from numpy import iterable,sum,array
        tree=self.getresp(self.makequery())
        if tree==None: return
        self.chanval[:,:]=0
        for ch in tree.findall("data"):
            if not "chan" in ch.attrib: continue
            id=int(ch.attrib["chan"])
            pixval=array([float(a) for a in ch.text.replace(",",".").split()])
            if id<len(self.chanval):
                self.chanval[id,:len(pixval)]=pixval
        return self.chanval
    
    def makesamp(self,data):
        '''creates Sample structure'''
        #from . import reband# import Sample
        from numpy import concatenate
        try:
            import reband
        except:
            from scanner import reband
        return reband.Sample(None,data=concatenate([[self.chanene[i],data[i]] for i in range(len(self.chanene))]))
    
    def setflat(self,minval=3,perc=3,fresh=True):
        from numpy import iterable,zeros_like
        if fresh: data=self.measure().copy()
        else: data=self.chanval.copy()
        if iterable(self.dark): data-=self.dark
        if (data<minval).sum()/data.size>0.5:
            print("intensity too low")
            return
        data[data<minval]=minval
        self.flat=data
        asel=zeros_like(data[0],dtype=bool)
        smooth=30
        asel[smooth:-smooth]=True
        #self.transtable,self.ysel=gettrans(self.pixtable,self.chanene,self.flat,smot=0.02,skiplowest=perc,osel=[asel]*3,disfun=lambda x:x**2,weifun=lambda x:x**3)
        self.transtable,self.ysel=gettrans(self.pixtable,self.chanene,self.flat,skiplowest=perc,
                                           osel=[asel]*3,disfun=lambda x:x**2,weifun=lambda x:x**3)

    def result(self,smooth=0,sub_dark=True,div_flat=True,maxit=0,scale=None):
        from numpy import sum,iterable
        if maxit==-2: data=self.chanval #reuse last measurement
        else: 
            rep=self.measure()
            if not iterable(rep):
                return None
            data=rep.copy()
        if sub_dark and iterable(self.dark): data-=self.dark
        if div_flat and iterable(self.flat):
                data/=self.flat
                if rc.hard_spec_min>0: data[data<rc.hard_spec_min]=rc.hard_spec_min
                if rc.hard_spec_max>0: data[data>rc.hard_spec_max]=rc.hard_spec_max
        if smooth>2:
            from numpy import hamming
            from scipy.ndimage.filters import convolve
            ham=hamming(smooth)
            ham/=ham.sum()
            if self.schan>=0:
                data=convolve(data,ham,mode="reflect")#mode="mirror")
            else:
                data=[convolve(dd,ham,mode="reflect") for dd in data]
        if iterable(scale) and len(scale)==len(data):
            for i in range(len(data)):
                data[i]*=scale[i]
        if maxit<0: return data
        if not iterable(self.ysel):
            print("to calibration first")
            return data
         #   self.transtable,self.ysel=gettrans(self.pixtable,self.chanene,self.flat)
        self.samp=self.makesamp(data)
        self.prof=self.samp.trans(self.transtable,self.ysel,lev=rc.band_minlev,ref=self.flat)
        self.last=self.prof.vals#sum([data[i][self.ysel[i]].dot(self.transtable[i]) for i in range(len(data))],axis=0)
        return self.last

class linkam():
    temprof=None
    def __init__(speed=-1,gotemp=None,ramp=30,init=True,close=False,port=0):
        import serial
        global ser
        if ser==None: ser=serial.Serial(port,baudrate=serspeed,timeout=2) # com1 communication
        if init:
            ser.write("\x3f\xa0\r")
            ser.readline()
        if gotemp: # not responding so far
            self.set_temp(speed,gotemp,ramp)
        if close:
            ser.close()
            ser=None
        #return temp
    def get_temp(self):
        ser.write("T\r")
        out=ser.readline()
        try:
            temp=int(out[-5:-1],16)*0.1
        except:
            print('error:'+out)
            return
        return temp
    def set_temp(self,speed=-1,gotemp=0,ramp=30):
        if speed>=0 and speed<=31:
            ser.write("P"+chr(48+speed)+"\r")
        ser.write("L1%03i\r"%int(gotemp*10))
        ser.write("R1%04i\r"%int(ramp*100))
        ser.write("Pa0\r")


    def reach_temp(self,temp,toler=1,per=2.):
        from time import sleep,time
        aa=linkam(gotemp=temp)
        at=aa.get_temp()
        ot=at
        t0=time()
        self.temprof=[]
        while abs(at-temp)>toler:
            lt=at
            sleep(per)
            aa.set_temp(gotemp=temp)
            at=aa.get_temp()
            if abs(lt-at)<toler/5.:
                print("not working")
                break
            t1=time()-t0
            if (lt-at)/per/(ot-at)*t1<0.1:
                print("cooling too slow")
                break
            self.temprof.append((t1,at))
        return at
