from traits.api import *
from traitsui.api import *
from threading import Thread
from time import sleep
from traitsui.key_bindings import KeyBinding, KeyBindings
from traitsui.qt4.extra.bounds_editor import BoundsEditor
from traitsui.ui_editors.array_view_editor import ArrayViewEditor

from traitsui.menu import Menu, Action, Separator
from traitsui.message import message
#from mpl_figure_editor import MPLFigureEditor
import matplotlib
#from scipy import indices,rand
from numpy import exp,sqrt,sum
global aw
aw=None

'''
import wx
matplotlib.use('wxAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wx import NavigationToolbar2Wx

from traitsui.wx.editor import Editor
from traitsui.wx.basic_editor_factory import BasicEditorFactory

'''

boost=1#50
from pyface.qt import QtGui, QtCore
from pyface.progress_dialog import ProgressDialog as ProgressBar

matplotlib.use('QT4Agg')
from matplotlib.figure import Figure

from scanner import labin,labin3d
from scanner import labrc as rc



'''
Event acquire - control pannel plots point in the design tab
acquisition thread : prepare, acquire, stack/...
problem with repeated callibration
'''

from scanner.my_editors import MPLFigureEditor, MPLInitHandler

#----------------------------------

class Experiment(HasTraits):
    """ Object that contains the parameters that control the experiment,
    modified by the user.
    """
    config=None
    instr=None
    stack={}
    reflist = ("unity","silicon","silicon+oxide (VASE)")

    expo = Float(rc.intime, label="exposure", desc="integ. time in miliseconds")
    aver = Int(rc.mean, label="averaging",desc="repeated meas.")
    smooth = Range(1,50,rc.smoothing,label="smooth",desc="level of smoothing")
    median = Int(1,label="software median",desc="more robust averaging")
    shut = Bool(True,label="light shutter",desc="open when checked")
    record = Bool(False,label="store next",desc="put in stack")
    combine = Bool(rc.chan_combine,label="combine multichannel")
    sname = File(label="Filename")
    refer = Button("Reference")
    darken = Button("Dark")
    saveme = Button("Save")
    saveall = Button("Save All")
    refermat = Enum(reflist)
    ready = Bool(False)
    paren = None
    errb = Bool(False,label="error",desc="show error band")

    menubar = MenuBar(Menu(Action(name='test Action', action='_run_action'), name="Menu"))

    view = View(Item('expo',width=5), Item('aver',resizable=True,width=5),
                HGroup(Item('shut'),Item('median',width=5),Item('combine'), Item('smooth',enabled_when="combine==False")),
                HGroup(Item('refer',show_label=False, enabled_when="ready"),Item('darken',show_label=False, enabled_when="ready")),
                HGroup(Item('sname'),Item('saveme',show_label=False, enabled_when="sname!=''"),Item('record')),
                HGroup(Item('saveall',show_label=False, enabled_when="sname!=''"),Item('refermat',label="Ref. material")),Item('errb'),#editor=CheckListEditor(values=reflist)
                menubar=menubar,width=10)

    def _shut_changed(self):
        import urllib2
        if self.instr==None: return
        query=self.instr.path+"0&shut=%i"%int(self.shut)
        #rep=urllib2.urlopen(query).read()
        message(query)

    def _expo_changed(self):
        if self.instr!=None: self.instr.intime=int(self.expo)

    def _combine_changed(self):
        if self.combine: self.smooth=1
        else: self.smooth=rc.smoothing

    def _aver_changed(self):
        if self.instr!=None: self.instr.config.m_NrAverages=int(self.aver)

    def _errb_changed(self):
        if not hasattr(self.instr,"prof"): return
        if self.errb:
            bmin,bmax=self.instr.prof.vals-2*self.instr.prof.errs,self.instr.prof.vals+2*self.instr.prof.errs
            self.bas=self.paren.figure.axes[0].fill_between(self.paren.spectrac.pixels, bmin,bmax,alpha=0.2)[0]
        else:
            if hasattr(self,'bas'):
                self.bas.remove()
                del self.bas

    def stack_name(self,defname='last'):
        if self.paren.scanner!=None and self.paren.scanner.ready:
            sname="pos"
            if len(self.paren.scanner.program)>0: sname="scan"
            sname+="_%03i_%03i"%tuple(self.paren.scanner.actpos)
            return sname
        if not self.record: return defname
        lstack=[i for i in self.stack.keys() if i.find(defname+"_")==0]
        self.record=False
        return defname+"_%05i"%(len(lstack)+1)

    def show_stack(self):
        print(" ".join(self.stack.keys()))

    def clear_stack(self):
        print(" ".join(self.stack.keys()))
        self.stack.clear()
        self.paren.stack_list.clear()

    def _saveme_fired(self):
        self.display(self.sname+" saved ...")
        from numpy import savetxt,array
        if self.combine==False: # all channels separate, not for stacks
            data=self.instr.result(maxit=-2)
            olist=[]
            for i in range(len(data)):
                olist+=[self.instr.chanene[i],data[i]]
            odata=array(olist)
        elif len(self.stack)>1:
            olist=self.stack.keys()
            odata=array([self.instr.pixtable]+[self.stack[k] for k in olist]).T
            print("saved %i spectra"%len(olist))
        else:
            odata=array([self.instr.pixtable,self.instr.last]).T
        savetxt(self.sname,odata,fmt="%8.5f")
        print(self.sname+" saved "+str(odata.shape))

    def _saveall_fired(self):
        from numpy import concatenate,array,savetxt
        self.display("stack saved ...")
        speclst=[concatenate([[0,0],self.instr.pixtable])]
        for k in self.stack.keys():
            scode=k.split("_")
            try:
                ix,iy=float(scode[1]),float(scode[2])
            except:
                continue
            #poslst.append([ix,iy])
            speclst.append(concatenate([[ix,iy],self.stack[k]]))
        savetxt(self.sname,array(speclst).T,fmt="%8.5f")
        print(self.sname+" saved /"+str(len(speclst)))

    def _darken_fired(self):
        if self.paren.acquisition_thread and self.paren.acquisition_thread.isAlive():
            self.paren.acquisition_thread.wants_abort = True
        message(" block the light beam, please ", title = 'User request')
        self.instr.dark=self.instr.result(sub_dark=False,div_flat=False,smooth=self.smooth,maxit=-1) #no correction/calibration
        self.paren.status_string="dark curr. saved"
        #if self.config!=None: self.config.adjust_image()

    def _refer_fired(self):
        if self.paren.acquisition_thread and self.paren.acquisition_thread.isAlive():
            self.paren.acquisition_thread.wants_abort = True
        self.paren.status_string="measuring reference sample"
        from numpy import any,iterable
        #if hasattr(self.instr,'ysel'): del self.instr.ysel
        rep=self.instr.result(div_flat=False,smooth=self.smooth,maxit=-1) #should use dark subtraction
        if not iterable(rep):
            print("reference failed")
            return
        self.instr.flat=rep
        if rc.saturate>0:
            for i in range(len(self.instr.flat)):
                if any(self.instr.flat[i]>=rc.saturate):
                    print("channel %i saturated"%i)
        if hasattr(self.instr,"setflat"):
            self.instr.setflat(perc=rc.flat_above_quantile,fresh=False)
            self.display("reference table calculated")
            self.display("channels :"+' '.join(["[%i]"%sum(a) for a in self.instr.ysel]))
            self.display("weights :"+' '.join(["[%.2f / %i]"%(a.sum(),sum(a.sum(0)>0)) for a in self.instr.transtable]))

        else: self.display("reference stored")
        self.paren.status_string="calibration finished"
        if not self.combine:
            for i in range(len(self.instr.flat)):
                self.instr.ysel[i]*=self.instr.flat[i]>=rc.underate
        if self.config!=None: self.config.adjust_image()

intwid=50
class Scan(HasTraits):

    instr=None

    Xpts = Int(10, label="X points", desc="fast axis/radial")
    Ypts = Int(10, label="Y points", desc="slow axis/azimuth")
    Xstep = Int(rc.cart_step, label="X step", desc="instrumental steps")
    Ystep = Int(rc.cart_step, label="Y step", desc="instrumental steps")
    centerstart = Bool(True, label="start at center")
    retlast = Bool(False, label="return to origin")
    zigzag = Bool(True, label="scan in both directions")
    radius = Int(50, label="wafer radius", desc="in mm")
    design = Button("Plan")
    dump = Button("Show")
    pclear = Button("Clear")
    setcenter = Button("Set center")
    cmeasure = Button("Find center")
    npoints=Int(0,label="Planned points ", enabled_when="program")
    speed=Range(1,10,label='Speed')
    cpos=String("[0,0]",label="Position")
    centstr=String("[0,0]",label="Center")
    centpos=[0,0]
    ready=Bool(False)
    gpos = Button("Get Position")

    up = Button("Up")
    down = Button("Down")
    left = Button("Left")
    right = Button("Right")
    shome = Button("Reset axis")
    gocenter = Button("Go center")
    srefer = Button("Calib sample")
    ptskip = Button("Skip")
    #res=[Button("%i s"%s) for s in [2,5,10]]
    view = View(VFold(
                HGroup(Item('Xpts',width=intwid),Item('Ypts',width=intwid),Item('Xstep',width=intwid),Item('Ystep',width=intwid)),
                #HGroup(Item('Xstep'),Item('Ystep')),
                HGroup(Item('radius', enabled_when="centerstart"),Item('centstr',style='readonly')),
                HGroup(Item('centerstart'),Item('retlast'),Item('zigzag'),Item('design',show_label=False),Item('dump',show_label=False),Item('pclear',show_label=False)),
                HGroup(Item('up',show_label=False),Item('down',show_label=False),
                    Item('left',show_label=False),Item('right',show_label=False),Item('ptskip',show_label=False),Item('gpos',show_label=False),enabled_when="ready"),
                HGroup(Item('npoints',style='readonly'),Item('cpos',style='readonly'),Item('shome',show_label=False),Item('srefer',show_label=False),Item('gocenter',show_label=False),Item('setcenter',show_label=False),
                    Item('cmeasure',show_label=False),enabled_when="ready"),
                Item('speed',show_label=True),
                label="Scanner", layout='split',
                )
            )
    program = []
    actpos = List(Float,[0,0])
    count = [0,0]
    exelist = {}

    def setup(self):
        if self.instr==None:
            self.instr=labin.specscope()
            self.instr.adrive(rc.motor_lin_rate,rc.motor_lin_gap)# may not be present
        self.ready=True
        self.centpos=rc.xy_cent
        self.centstr=str(list(self.centpos))

    def _design_fired(self):
        from numpy import mgrid,newaxis,array
        self.centpos=rc.xy_cent#self.setup()
        self.centstr=str(list(self.centpos))
        self.program=mgrid[:self.Xpts,:self.Ypts]*(array([self.Xstep,self.Ystep]).reshape(2,1,1))
        if self.zigzag:
            for i in range(1,self.Xpts,2):
                self.program[1][i]=self.program[1][i][::-1]
        center=array(self.centpos).reshape(2,1,1)
        if self.centerstart:
            self.program+=center
            self.program[0]-=int(self.Xpts*self.Xstep/2.)
            self.program[1]-=int(self.Ypts*self.Ystep/2.)
        if self.radius>0:
            #for i in range(2):
            #    self.program[i]+=self.centpos[i]
            #if self.radius>0:
            #self.program=array(self.program)
            sel=((self.program-center)**2).sum(0)<self.radius**2
            print("selected %i points"%sum(sel))
            self.program=list(self.program[:,sel].T)
        else:
            self.program=list(self.program.reshape(2,self.Xpts*self.Ypts).T)
        if self.retlast:
            self.program.append(self.program[0])
        self.npoints=len(self.program)
        if 'plan' in self.exelist: #vykreslit
            self.exelist['plan'].design_show(self.program)
            self.exelist['plan'].experiment.clear_stack()

    def _srefer_fired(self):
        self.instr.goto(rc.refer_pos[0],rc.refer_pos[1])

    def _shome_fired(self):
        if self.instr: self.instr.ahome()
        for i in range(2):
            self.actpos[i]=rc.xy_cent[i]
        self.cpos=str(list(self.actpos))
        self.setup()

    def _gpos_fired(self):
        self.instr.awrite("M114")
        inp=self.instr.acomm()
        for l in inp:
            sl=str(l)
            if sl.find('X')>0:
                #print("got resp "+sl)
                for k in sl.split('.0'):
                    if k.find('X')>=0:
                        ival=float(k[k.find(':')+1:])
                        self.actpos[0]=int(ival)
                    elif k.find('Y')>=0:
                        ival=float(k[k.find(':')+1:])
                        self.actpos[1]=int(ival)
                self.cpos=str(list(self.actpos))

    def _cmeasure_fired(self):
        print("trying to find wafer center [currently %s]"%str(self.centpos))
        self.centpos=self.actpos
        self.centstr=str(list(self.centpos))

    def _setcenter_fired(self):
        print("current position as wafer center")
        self.centpos=self.actpos
        self.centstr=str(list(self.centpos))

    def _actpos_changed(self):
        self.cpos=str(list(self.actpos))

    def _gocenter_fired(self):
        self.instr.goto(self.centpos[0],self.centpos[1])

    def _right_fired(self):
        if self.instr: self.instr.rate(2,self.Ystep)
        self.actpos[1]+=self.Ystep

    def _left_fired(self):
        if self.instr: self.instr.rate(2,-self.Ystep)
        self.actpos[1]-=self.Ystep

    def _down_fired(self):
        if self.instr: self.instr.rate(1,-self.Xstep)
        self.actpos[0]-=self.Xstep

    def _up_fired(self):
        if self.instr: self.instr.rate(1,self.Xstep)
        self.actpos[0]+=self.Xstep

    def _ptskip_fired(self):
        if len(self.program)==0: return
        self.actpos=self.program.pop(0)
        self.cpos=str(list(self.actpos))
        if 'plan' in self.exelist:
            self.exelist['plan'].acquire=True

    def _dump_fired(self):
        if len(self.program)>0:
            from numpy import array
            print(array(self.program))
        print("position: %i %i"%tuple(self.actpos))

    def _pclear_fired(self):
        self.program=[]
        if 'plan' in self.exelist: #vykreslit
            self.exelist['plan'].design_show(self.program)
            self.exelist['plan'].experiment.clear_stack()
        self.npoints=len(self.program)

    def next_point(self):
        from time import sleep
        if len(self.program)==0: return
        point=self.program.pop(0)
        if hasattr(self.instr,'ard'):
            if hasattr(self.instr,'goto'):
                self.instr.goto(point[0],point[1])
            else:
                self.instr.rate(1,point[0]-self.actpos[0])
                self.instr.rate(2,point[1]-self.actpos[1])
        else:
            print("riding %i,%i"%(point[0],point[1]))
            sleep(10-self.speed)
        self.actpos=list(point)
        self.npoints-=1
        #self.cpos=str(list(self.actpos))
        for k in self.exelist.keys():
            self.exelist[k].acquire=True

    def run_scan(self):
        #obsolete
        #now part of acquisition thread
        self.process=self.instr.scan_save(self.Xstep,self.Ystep,self.Xpts,self.Ypts,self.radius)
#----------------------------------------

def analyse(spect):
    """ Function called to do the processing (more time consuming, doesn't block the acquisition) - nothing implemented yet"""
    from scanner import spectra
    #spectra.extrema()

class AcquisitionThread(Thread):
    """ Acquisition loop. This is the worker thread that retrieves data, displays them, and spawns the processing job.
    """
    wants_abort = False
    n_img = 0

    def process(self, spect):
        """ Spawns the processing(analysis) job. """
        try:
            if self.processing_job.isAlive():
                self.display("Processing too slow")
                return
        except AttributeError:
            pass
        self.analyse.calculate = True
        self.processing_job = Thread(target=analyse, args=(spect,))#,self.results))
        self.processing_job.start()
        self.result()

    def run(self):
        """ Runs the acquisition loop. """
        self.display('Spectrac started')
        self.n_img = 0
        from numpy import iterable
        while not self.wants_abort:
            self.n_img += 1
            if hasattr(self,"prepare"): self.prepare()
            spect=self.acquire(self.experiment)
            if not iterable(spect):
                print("acquisition failed")
                break #returned none
            if self.experiment.refermat!='unity':
                instr=self.experiment.instr
                if instr.samp!=None:
                   instr.samp.calib()
                   instr.samp.prof=instr.samp.trans(instr.transtable,instr.ysel,lev=rc.band_minlev,ref=instr.flat)
                   spect=instr.samp.prof.vals
                else:
                    print("sample not created")
            self.display('%d spectra captured' % self.n_img)
            if hasattr(self,"stack"):
                sname=self.experiment.stack_name()
                if not sname in self.stack:
                    self.experiment.paren.stack_list.append(sname)
                self.stack[sname]=spect
            self.image_show(spect)
            self.process(spect)
        self.display('Spectrac stopped')

#--------------------------------------------------------------------
class Spectrac(HasTraits):
    """ Spectrac setup and acquisition.
    """
    erange = Range(0.,1.0)
    #erange = grange
    elow = Float(rc.meas_range[0], label="from", desc="lower energy band")
    ehigh = Float(rc.meas_range[1], label="to", desc="upper energy band")
    #gain = Enum(1, 2, 3, label="Gain", desc="gain")
    simulate = Bool(False, label="Emulate")
    port = Int(rc.web_port, label="connect. port no.")
    #Item('elow'), Item('ehigh'),
    setmeup = Button("Setup")
    equalize = Button("Equal")
    calib = Button("Calibrate Si")
    chanshow = Button("Show channels")
    resol = Float(rc.meas_comb_resol, label="spectral resolution in eV")
    cooler=Bool(False)
    cooltemp = Range(-20,10,-10,label="Cooling",desc="required temp. in centigrade")
    norm = Tuple(Float,Float,Float)#List(Float, [1.,1,1])
    darkcorr = Bool(True, label="Dynamic dark")
    nonlincorr = Bool(True, label="Nonlin. correction")

    chanmatch = Button("Match")
    refer_x = Int(rc.refer_pos[0])
    refer_y = Int(rc.refer_pos[1])
    #norm = Array
    view = View(VGroup(
                Item('erange', editor=BoundsEditor(low_name = 'elow', high_name = 'ehigh')),
                Item('simulate'),
                HSplit(Item('resol',width=5),Item('cooltemp',width=5, enabled_when='cooler'),),
                HSplit(Item('darkcorr',width=5),Item('nonlincorr',width=5),),
                HSplit(Item('norm',style='simple',label="channel equalizer", enabled_when='instr!=None'),Item('equalize',width=5,show_label=False),
                    Item('chanshow',width=5,show_label=False),Item('chanmatch')),#,enabled_when="instr.flat!=None")),
                Item('port'),Item('calib',width=5,show_label=False, enabled_when='instr.samp!=None'),
                HSplit(Item('refer_x'),Item('refer_y'),),
                Item('setmeup',show_label=False),
                springy=True),resizable=True)#, enabled_when="config!=None")
    from numpy import array
    #norm = array([1,1.,1.])
    cooler=False
    instr=None
    pixels=None
    exper=None
    paren=None

    acquire = Event
    def _chanmatch_fired(self):
        if not hasattr(self.instr,'samp') or self.instr.samp==None: return
        corrs=self.instr.samp.mismat(self.instr.transtable)
        print(corrs)
        for i,b in enumerate(self.instr.samp.bands):
            if i>0:
                b.scale=1/corrs[i-1].n

    def setup(self):
        #for i in range(3):
        #    self.norm[i]=1.
        if self.instr==None:
            if self.simulate: self.instr=labin.specscope()
            elif hasattr(rc,"xy_cent"):
                self.instr=labin3d.webocean3d(path="http://localhost:%i/?exp="%self.port,chan=rc.chan_sel)
                self.instr.gxmax,self.instr.gymax=rc.xy_size
            else: self.instr=labin.uniocean(path="http://localhost:%i/?exp="%self.port,chan=rc.chan_sel) ##FM replaced webocean
        self.instr.setup([self.elow,self.ehigh],self.resol,integ=self.exper.expo,aver=self.exper.aver)
        if len(self.instr.chanene)>0: self.exper.display(" found %i channels ..."%len(self.instr.chanene))
        #if len(self.instr.chanrange)>0: self.exper.display(" found %i channels ..."%len(self.instr.chanrange))
        self.exper.instr=self.instr
        if hasattr(self.instr,'ard'):
            self.exper.display("scanner connected to %s"%rc.ard_port)
            self.paren.scanner.instr=self.instr
            self.paren.scanner.ready=True
        self.pixels=self.instr.pixtable.copy()#
        if self.instr.config.m_StartPixel>0 or (self.instr.config.m_StopPixel>1 and self.instr.config.m_StopPixel<len(self.pixels)):
            self.pixels=self.pixels[self.instr.config.m_StartPixel:self.instr.config.m_StopPixel]
        self.instr.config.Material=b'simu'
        self.instr.samp=None
        self.exper.display("setup ok [%i pixels] ..."%len(self.pixels))
        self.norm=(1,1,1)
        self.paren.status_string="now calibrate [Experiment/Reference+Dark]"

    def _setmeup_fired(self):
        self.setup()
        self.exper.ready=True
        self.paren.ready=True
        #self.scanner.setup()

    def _chanshow_fired(self):
        if self.instr==None: return
        #self.paren.figure.clean()
        self.paren.figure.axes[0].clear()
        self.paren.spect_last=[]
        for i in range(len(self.instr.chanene)):
           out=self.paren.figure.axes[0].plot(self.instr.chanene[i],self.instr.flat[i])[0]
           self.paren.spect_last.append(out)
        self.paren.figure.canvas.draw()

    def _norm_changed(self):
        self.instr.intfact=list(self.norm)

    def _equalize_fired(self):
        from numpy import iterable,array,percentile
        if not(iterable(self.instr.flat)): return
        if not(iterable(self.instr.flat[0])): return
        from numpy import array,round
        #zmax=array([c.max() for c in self.instr.flat])
        zmax=array([percentile(c,90) for c in self.instr.flat])
        self.exper.display("maxima [%s] ..."%str(zmax))
        zmax/=zmax.mean()
        zmax[zmax<=0]=1
        zmax=1./zmax
        if hasattr(self.instr,'transtable') and iterable(self.instr.transtable):
            tsum=[sum(a) for a in self.instr.transtable]
            for i in range(len(zmax)):
                if tsum[i]==0: zmax[i]=0
        #for i in range(len(zmax)):
        orig_norm=array(self.norm)
        zmax*=orig_norm #values used in previous acquisition
        self.norm=tuple(round(zmax,3))
        self.instr.intfact=list(self.norm)
        combchan=array([a.sum(0)>0 for a in self.instr.transtable]).sum(0)
        subpix=self.instr.pixtable[combchan>0]
        if len(subpix)>0:
            print("sum combined %i [%.2f - %.2f]"%(sum(combchan),subpix[0],subpix[-1]))
    #def _acquire_fired(self):
        #self.measure()
    #    self.config.display("arrived")

    def _calib_fired(self):
        print(self.exper.refermat)
        assert self.instr.samp!=None
        self.instr.samp.calib()
        print([b.samp.norm(b.ix).mean() for b in self.instr.samp.bands])

    def close(self):
        if self.instr!=None:
            self.instr.end()
            if rc.debug==0: del self.instr

    def measure(self,experiment=None): #main task during acquisition
        if self.instr==None:
            self.exper.display("must configure spectroscope first")
            return None
        self.instr.darkcorr=self.darkcorr
        self.instr.nonlincorr=self.nonlincorr
        if experiment.combine:
            if experiment.median>1:
                from numpy import median
                allspect=[self.instr.result() for i in range(experiment.median)]
                spect=median(allspect,0)
            else:
                spect=self.instr.result()
            if not(hasattr(spect,"shape")): return None
            if len(spect.shape)>1:
                return spect.sum(axis=0)
        else:
            spect=self.instr.result(maxit=-1,smooth=experiment.smooth)
        # smooth should be called internally as result(smooth=True)
        if (experiment!=None and experiment.combine and experiment.smooth>20):
            from numpy import hamming,convolve
            ham=hamming(experiment.smooth)
            ham/=ham.sum()
            spect=convolve(ham,spect,"same")
        if experiment.combine:
            from numpy import where
            selpix=where(spect>0)[0]
            if len(selpix>0):
                if selpix[0]>0: spect[:selpix[0]]=spect[selpix[0]]
                if selpix[-1]<len(spect)-1: spect[selpix[-1]+1:]=spect[selpix[-1]]
            else:
                print("problem combining:no pict selected")
        return(spect)

    def get_match(self):
        perc=rc.flat_above_quantile
        ou,ol=labin.gettrans(self.instr.pixtable,self.instr.chanene,self.instr.flat,skiplowest=perc,rep=-1)
        data=self.instr.result(maxit=-2) # reuse last measurement
        #WRONG: last data in chanval are not flat corrected
        #zval=[self.instr.chanval[i][ol[i]].dot(ou[i]) for i in range(len(ou))]
        zval=[data[i][ol[i]].dot(ou[i]) for i in range(len(ou))]
        istart,iend=overlaps(caltab)
        for i in range(len(istart)):
            cors=[zval[j][istart[i]:iend[i]].mean() for j in range(len(ou))]
            print("%.2f-%.2f:"%(istart[i],iend[i])+str(cors))

    def _refer_x_changed(self):
        rc.refer_pos=[int(self.refer_x),int(self.refer_y)]

    def _refer_y_changed(self):
        rc.refer_pos=[int(self.refer_x),int(self.refer_y)]

class Analyse(HasTraits):
    erange = Range(0.,1.0)
    #erange = grange
    elow = Float(0.9, label="from", desc="lower energy band")
    ehigh = Float(5.5, label="to", desc="upper energy band")
    code = Code("return xdat.mean(),ydat.mean()")
    run = Button('Eval')
    res = Button('Reset')
    output = String("here comes the data",label="Output")
    view = View(
                Item('erange', editor=BoundsEditor(low_name = 'elow', high_name = 'ehigh')),
                Item('code',show_label=False),
                HGroup(Item('run',show_label=False),Item('res',show_label=False),),
                Item('output',show_label=False,style='readonly')
                )#, enabled_when="config!=None")
    calculate = Event
    func=None
    vals=[]

    def evalme(self):
        xdat=self.instr.pixtable
        sel=(xdat>=self.elow)*(xdat<=self.ehigh)
        ydat=self.instr.last
        return self.func(xdat[sel],ydat[sel])

    def _reset_fired(self):
        self.vals=[]

    def _run_fired(self):
        #rep=exec(self.code) not working
        imp=self.code
        imp="def runme(xdat,ydat):\n    "+"\n    ".join(self.code.split("\n"))
        open("analme.py","w").write(imp)
        if 'analme' in locals():
            from imp import reload
            del self.func
            reload(analme)
        else: import analme
        self.func=analme.runme
        self.output=str(self.evalme())
        #self.output="Evaluation failed"

    def _calculate_fired(self):
        if self.func!=None:
            self.vals=self.evalme()
            self.output=str(self.vals)

my_bindings = KeyBindings(
    KeyBinding (binding1 = 'r',
                description = 'Start/Stop',
                method_name = '_start_stop_acquisition'),
        )

from traitsui.tabular_adapter import TabularAdapter
class MultiSelectAdapter(TabularAdapter):
    """ This adapter is used by both the left and right tables
    """
    columns = [('', 'myvalue')]
    myvalue_text = Property
    def _get_myvalue_text(self):
        return self.item

class ControlPanel(HasTraits):
    """ This object is the core of the traitsUI interface. Its view is
    the right panel of the application, and it hosts the method for
    interaction between the objects and the GUI.
    """
    experiment = Instance(Experiment, ())
    spectrac = Instance(Spectrac, ())
    scanner = Instance(Scan, ())
    figure = Instance(Figure)
    design = Instance(Figure)
    results = Instance(Figure)
    start_stop_acquisition = Button("Start/Stop")
    results_string = String()
    stack_list = List(Str, [])
    stack_selec = List(Str, [])
    analyse = Instance(Analyse, ())
    acquisition_thread = Instance(AcquisitionThread)
    #bardial=ProgressBar(min=0,max=100,title='acquired')
    #bar=bardial.progress_bar
    spect_last=None
    ready=Bool(False)
    #ready=False
    acquire = Event
    status_string = String();
    result_plot=[]
    grcolors="rgbkmc"
    showme=Button("Show")
    delme=Button("Delete")
    stack_disp=[]

    view = View(Group(
                    Group(
                        Group(
                            Item('experiment', style='custom', show_label=False),
                            Item('scanner', style='custom', show_label=False),
                            label="Input",),
                        label='Experiment', dock="tab"),
                    Group(
                        Item('results_string',show_label=False, springy=True, style='custom' ),
                        #Item('bar', show_label=False ),
                        UItem('stack_list', label='Stack',
                              editor=TabularEditor(show_titles=True,selected='stack_selec',editable=True,
                              multi_select=True,adapter=MultiSelectAdapter())
                            #editor=ListEditor(style='simple',rows=1)
                        ),
                        Group(Item("showme", show_label=False)),Group(Item("delme", show_label=False)),
                        label="Control", dock='tab',),
                    Group(
                            Item('spectrac', style='custom', show_label=False,),
                            #Separator(),

                        label='Machine',dock="tab", selected=True),
                    Group(
                        Group(
                            Item('analyse', style='custom', show_label=False),
                            label="Code",),
                        label='Analyse', dock="tab"),
                    layout='tabbed'),
                   # Group(
                        Item('status_string', show_label=False, style='readonly'),
                        Item('start_stop_acquisition', show_label=False, enabled_when="experiment.ready" ),

                   # ),
                resizable=True)

    def __init__(self,*args,**kwarks):
        HasTraits.__init__(self,*args,**kwarks)
        #self.results.display=self.figure
        self.experiment.display = self.add_line
        self.experiment.paren = self
        self.spectrac.exper = self.experiment
        self.spectrac.paren = self
        #if rc.auto_init: self.spectrac.setup()
        self.scanner.exelist['spectrac']=self.spectrac
        self.scanner.exelist['plan']=self
        #self.bardial.open()

    def _showme_fired(self):
        from matplotlib import patches
        for p in self.stack_disp:
            p.remove()
            del p
        self.stack_disp=[]
        i=0
        for a in self.stack_selec:
            if not(a) in self.experiment.stack: continue
            #print(a,sum(self.experiment.stack[a]))
            if self.experiment.combine:
                mcol=rc.all_color[i%len(rc.all_color)]
                p=self.figure.axes[0].plot(self.spectrac.pixels,self.experiment.stack[a],mcol)[0]
                self.stack_disp.append(p)
                coors=a.split("_")
                i+=1
                if len(coors)==3:
                    try:
                        qx,qy=int(coors[1]),int(coors[2])
                        p=patches.CirclePolygon(qx,qy,1,color=mcol)
                        self.stack_disp.append(p)
                    except:
                        print("failed parse of "+a)
                        pass
        #print("%i items shown"%(len(self.stack_disp)))
        self.figure.canvas.draw()
        self.design.canvas.draw()

    def _delme_fired(self):
        for a in self.stack_selec:
            i=self.stack_list.index(a)
            if i>=0:
                del self.stack_list[i]
                if a in self.experiment.stack:
                    del self.experiment.stack[a]


    def _start_stop_acquisition_fired(self):
        """ Callback of the "start stop acquisition" button. This starts
        the acquisition thread, or kills it.
        """
        if self.acquisition_thread and self.acquisition_thread.isAlive():
            self.acquisition_thread.wants_abort = True
            for k in self.experiment.stack.keys():
                if not k in self.stack_list:
                    self.stack_list.append(k)
            self.status_string="Stopped"
            #self.start_stop_acquisition.label_value="Stopped"
        else: #starting....
            #self.figure.clean()
            from numpy import iterable
            if self.spect_last!=None:
                if iterable(self.spect_last):
                    cnt=len(self.spect_last)
                    for sp in self.spect_last: sp.remove()
                else:
                    self.spect_last.remove()
                self.spect_last=None
            if self.spectrac.instr==None:
                message("press Setup first", title = 'User request')
                return
            self.scanner.instr=self.spectrac.instr
            self.analyse.instr=self.spectrac.instr

            self.acquisition_thread = AcquisitionThread()
            if self.scanner.ready: #scanning in process
                self.acquisition_thread.prepare = self.scanner.next_point
            self.acquisition_thread.display = self.add_line
            self.acquisition_thread.acquire = self.spectrac.measure
            self.acquisition_thread.experiment = self.experiment
            self.acquisition_thread.stack = self.experiment.stack
            self.acquisition_thread.image_show = self.image_show
            self.acquisition_thread.analyse = self.analyse
            self.acquisition_thread.result = self.result_update
            #self.bardial.update(10)
            #self.start_stop_acquisition.label_value="Running"
            self.status_string="Running"
            self.acquisition_thread.start()

    def adjust_image(self,margin=0):
        from numpy import iterable
        if self.spect_last!=None:
            if iterable(self.spect_last):
                glist=self.spect_last
            else:
                #xval=[self.spectrac.instr.pixtable]
                glist=[self.spect_last]
            rang=self.figure.axes[0].get_xlim()
            glist=self.figure.axes[0].lines #pokus
            ymin,ymax=1000.,-10.
            for i in range(len(glist)):
                xval=glist[i].get_xdata()
                sel=(xval>=rang[0])*(xval<=rang[1])
                if sum(sel)<3: continue
                yval=glist[i].get_ydata()[sel]
                if ymin>yval.min(): ymin=yval.min()
                if ymax<yval.max(): ymax=yval.max()
            if margin>0:
                ydis=ymax-ymin
                ymin,ymax=ymin-ydis*margin,ymax+ydis*margin
            if ymin<rc.graph_min: ymin=rc.graph_min
            if ymax>rc.graph_max: ymax=rc.graph_max
            self.figure.axes[0].set_ylim(ymin,ymax)
            self.figure.canvas.draw()

    def add_line(self, string):
        """ Adds a line to the textbox display."""
        self.results_string = (string + "\n" + self.results_string)[0:1000]

    def image_show(self, spect):
        """ Plots an image on the canvas in a thread safe way. """
        #self.figure.axes[0].images=[]
        if self.spect_last==None:
            if rc.debug>0: print("(re)new graphs")
            if self.experiment.combine:
                self.spect_last=self.figure.axes[0].plot(self.spectrac.pixels,spect,rc.line_color,lw=rc.line_width)[0]
            else:
                pin=self.spectrac.instr
                self.spect_last=[]
                for i in range(len(spect)):
                    dx,dy=pin.chanene[i],spect[i]
                    if pin.ysel!=None:
                        dx,dy=dx[pin.ysel[i]],dy[pin.ysel[i]]
                    self.spect_last.append(self.figure.axes[0].plot(dx,dy,lw=rc.line_width,color=self.grcolors[i])[0])
            self.figure.axes[0].grid(1)
        else:
            if self.experiment.combine:
                self.spect_last.set_ydata(spect)
            else:
                pin=self.spectrac.instr
                for i in range(len(spect)):
                    if self.spectrac.norm[i]>0: #instr.intfact
                        self.spect_last[i].set_ydata(spect[i][pin.ysel[i]])
        self.figure.canvas.draw()
        #self.profile.canvas.draw()

    def result_update(self):
        """plots results of analysis"""
        from numpy import arange,iterable
        if iterable(self.analyse.vals):
            vals=self.analyse.vals
            nvals=len(vals)
        else:
            vals=[self.analyse.vals]
            nvals=1
        if nvals==0: return
        if len(self.results.axes)==0:
            divs=arange(0.02,0.99,0.96/nvals)
            for i in range(nvals):
                self.results.add_axes([0.05,divs[i],0.92,0.96/nvals-0.06])
            for i in range(nvals):
                ax=self.results.axes[i]
                self.result_plot.append(ax.plot(1,vals[i])[0])
                ax.grid()
                x=[1]
        else:
            for i in range(len(self.results.axes)):
                pl=self.result_plot[i]
                x,y=pl.get_xdata(),pl.get_ydata()
                from numpy import concatenate
                x=concatenate([x,[len(x)]])
                y=concatenate([y,[vals[i]]])
                pl.set_xdata(x)
                pl.set_ydata(y)

                ax=self.results.axes[i]
                ax.set_xlim(max(0,len(x)-40),len(x)+1)
                yran=[y.min(),y.max()]
                ax.set_ylim(yran[0]*1.1-yran[1]*0.1,yran[1]*1.1-yran[0]*0.1)
        self.results.canvas.draw()

    def design_show(self, plan):
        ax=self.design.axes[0]
        ax.clear()
        if len(plan)>0:
            ax.plot([p[0] for p in plan],[p[1] for p in plan],'bd-')
            pp=ax.get_xlim()
            mdis=(pp[1]-pp[0])*0.02
            ax.set_xlim(pp[0]-mdis,pp[1]+mdis)
            pp=ax.get_ylim()
            ax.set_ylim(pp[0]-mdis,pp[1]+mdis)
        self.design.canvas.draw()

    def _acquire_fired(self):
        self.design.axes[0].plot(self.scanner.actpos[0],self.scanner.actpos[1],'rs')
        self.design.canvas.draw()

class MainWindowHandler(Handler):
    def close(self, info, is_OK):
        #sleep(1)
        if ( info.object.panel.acquisition_thread
            and info.object.panel.acquisition_thread.isAlive() ):
            info.object.panel.acquisition_thread.wants_abort = True
            while info.object.panel.acquisition_thread.isAlive():
                sleep(0.1/boost)
        #    taking care of instruments
        info.object.panel.spectrac.close()
        return True

class MainWindow(HasTraits):
    """ The main window, here go the instructions to create and destroy the application. """
    #display = Instance(Display)

    figure=Instance(Figure)
    design=Instance(Figure)
    results=Instance(Figure)
    panel = Instance(ControlPanel)

    adjust_image = Button("Adjust graph")

    def _figure_default(self):
        figure = Figure()
        figure.add_axes([0.1, 0.04, 0.87, 0.92])
        return figure

    def _design_default(self):
        design = Figure()
        design.add_axes([0.1, 0.04, .87, .92])
        design.axes[0].grid()
        return design

    def _results_default(self):
        results = Figure()
        #results.add_axes([0.05, 0.04, .9, .92])
        #results.axes[0].grid()
        return results

    def _adjust_image_fired(suself):
        suself.panel.adjust_image()

    def _panel_default(self):
        return ControlPanel(figure=self.figure,design=self.design,results=self.results)

tabbedview = View(HSplit(Group(Group(
                                Item('figure', editor=MPLFigureEditor(), label="Spectrum", dock='tab',resizable=True),
                                Item('design', editor=MPLFigureEditor(), dock='tab'),
                                Item('results', editor=MPLFigureEditor(), dock='tab'),
                                layout="tabbed", show_labels=False,springy=True),
                            Item('adjust_image', show_label=False ),
                            springy=True),
                       Item('panel', style="custom",resizable=True),
                       show_labels=False,
                      ),
                resizable=True,
                height=0.75, width=0.85,
                handler=MainWindowHandler(),
                title='Spectrac NG',
                key_bindings = my_bindings,
                buttons=NoButtons)


#--------------------extra stuff-----------------
def fm_init(self):
    import pyfirmata as fm
    self.board = pyfirmata.Arduino(self.port)

    #do we need this? probably only for inputs
    it = pyfirmata.util.Iterator(self.board)
    it.start()

    self.board.digital[int(self.pin)].write(0)

    value=float(self.board.analog[int(self.pin2)].read())

    digital_0 = self.board.get_pin('d:11:p') #pwm
    digital_0.write(value) #between 0 and 1


import sys
if __name__ == '__main__':
    aw=MainWindow()
    #if "--init" in sys.argv: rc.auto_init=True
    aw.configure_traits(view=tabbedview)
    message(" press Setup first ", title = 'User request')
        