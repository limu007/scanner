#from traits.api import *
from traitsui.api import View, Group, HGroup, HSplit, VGroup, VFold, NoButtons
from traitsui.api import EnumEditor, TabularEditor, FileEditor, Menu, MenuBar, Item, UItem, Handler
from traits.api import HasTraits, Int, Float, Bool, Str, Range, String, List, Tuple, Enum, File
from traits.api import Event, Code, Button, Property, Instance

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
from numpy import exp,sqrt,sum, array, iterable
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
from matplotlib import pyplot#.cm import get_cmap
glob_cmap=pyplot.cm.get_cmap('RdYlBu')

from . import labin,labin3d
from . import labrc as rc



'''
Event acquire - control pannel plots point in the design tab
acquisition thread : prepare, acquire, stack/...
problem with repeated callibration
'''

from scanner.my_editors import MPLFigureEditor, MPLInitHandler, SelectFromCollection

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
    shut = Bool(rc.use_shut,label="light shutter",desc="open shutter for measurement only")
    record = Bool(False,label="store next",desc="put in stack")
    combine = Bool(rc.chan_combine,label="combine multichannel")
    savecols = Bool(True,label="spectra in columns")
    sname = File(label="Filename")
    refer = Button("Reference")
    darken = Button("Dark")
    saveme = Button("Save")
    saveall = Button("Save All")
    refermat = Enum(reflist)
    ready = Bool(False)
    paren = None
    errb = Bool(False,label="errorband",desc="whether display error band")
    recalper = Int(0,label="recalib. period",desc="regular calibration during scanning")
    material = String("SiO2/Si",label="layer. structure")

    menubar = MenuBar(Menu(Action(name='test Action', action='_run_action'), name="Menu"))

    view = View(HGroup(Item('expo',width=5), Item('aver',resizable=True,width=5),),
                HGroup(Item('shut'),Item('median',width=5),Item('combine'), Item('smooth',enabled_when="combine==False")),
                HGroup(Item('darken',show_label=False, enabled_when="ready"),Item('refer',show_label=False, enabled_when="ready"),
                Item('recalper',enabled_when="ready"),),
                HGroup(Item('sname'),Item('saveme',show_label=False, enabled_when="sname!=''"), Item('record')),
                HGroup(Item('saveall',show_label=False, enabled_when="sname!=''"),Item('savecols'), Item('refermat',label="Ref. material")),
                Item('errb'),#editor=CheckListEditor(values=reflist)
                menubar=menubar,width=10)

    def _shut_changed(self):
        #import urllib2
        if self.instr==None: return
        rc.use_shut=1 if self.shut else 0
        #query=self.instr.path+"0&shut=%i"%int(self.shut)
        #rep=urllib2.urlopen(query).read()
        #message(query)

    def _expo_changed(self):
        if self.instr!=None: 
            self.instr.intime=int(self.expo)
            self.instr.set_acquisition(self.expo,self.aver)

    def _combine_changed(self):
        if self.combine: self.smooth=1
        else: self.smooth=rc.smoothing

    def _aver_changed(self):
        if self.instr!=None: 
            self.instr.set_acquisition(self.expo,self.aver)
            #self.instr.config.m_NrAverages=int(self.aver)

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
            sname+="_%03i_%03i"%tuple(self.paren.scanner.actpos[:2])
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
                if hasattr(self,'intfact') and self.intfact[i]==0: continue
                olist+=[self.instr.chanene[i],data[i]]
            odata=array(olist).T
        elif len(self.stack)>1:
            olist=self.stack.keys()
            odata=array([self.instr.pixtable]+[self.stack[k] for k in olist])
            print("saved %i spectra"%len(olist))
        else:
            odata=array([self.instr.pixtable,self.instr.last])
        if not self.savecols: odata=odata.T
        #todo check for existing filename
        #incremental naming
        savetxt(self.sname,odata,fmt="%8.5f")
        print(self.sname+" saved "+str(odata.shape))

    def _saveall_fired(self,selected=False):
        from numpy import concatenate,array,savetxt
        self.display("stack saved ...")
        if self.combine: clist=[0]
        else: 
            clist=[i+1 for i in range(len(self.instr.chanene)) if self.instr.intfact[i]>0]
        for cn in clist: #looping through channel numbers
            if cn>0: speclst=[concatenate([[0,0],self.instr.chanene[cn-1]])]
            else: speclst=[concatenate([[0,0],self.instr.pixtable])]
            for k in self.stack.keys():
                if selected and not k in self.paren.stack_selec:
                    continue
                scode=k.split("_")
                try:
                    ix,iy=float(scode[1]),float(scode[2])
                except:
                    continue
                #poslst.append([ix,iy])
                if cn>0 and (type(self.stack[k])==list or len(self.stack[k].shape)>1):
                    speclst.append(concatenate([[ix,iy],self.stack[k][cn-1]]))
                else:
                    speclst.append(concatenate([[ix,iy],self.stack[k]]))
            oname=self.sname
            if cn>0:
                if oname.find('.')>0: oname=oname.replace(".","_%i."%cn)
                else: oname=oname+"_%i."%cn
            try:
                savetxt(oname,array(speclst).T,fmt="%8.5f")
            except:
                print("cannot save")
            print(oname+" saved /"+str(array(speclst).shape))

    def _darken_fired(self,interrupt=True):
        from numpy import median, iterable
        if iterable(self.instr.dark): odark=median(self.instr.dark,1)
        else: odark=[]
        if self.paren.acquisition_thread and interrupt and self.paren.acquisition_thread.isAlive():
            self.paren.acquisition_thread.wants_abort = True
        if self.shut:
            rc.use_shut=0
        else:
            message(" block the light beam, please ", title = 'User request')
        self.paren.status_string="dark acquisition"
        self.instr.dark=self.instr.result(sub_dark=False,div_flat=False,smooth=self.smooth,maxit=-1) #no correction/calibration
        self.paren.status_string="dark curr. saved"
        if len(odark)==len(self.instr.dark):
            odark=odark-median(self.instr.dark,1)
            self.display("dark measured - dif %s"%str(odark))
        else:
            odark=median(self.instr.dark,1)
            self.display("dark measured - median %s"%str(odark))
        if self.shut:
            rc.use_shut=1
        #if self.config!=None: self.config.adjust_image()

    def _refer_fired(self,interrupt=True):
        if self.paren.acquisition_thread and interrupt and self.paren.acquisition_thread.isAlive():
            self.paren.acquisition_thread.wants_abort = True
        self.paren.status_string="calibration running"
        self.paren.status_string="measuring reference sample"
        from numpy import any,iterable,percentile
        #if hasattr(self.instr,'ysel'): del self.instr.ysel
        rep=self.instr.result(div_flat=False,smooth=self.smooth,maxit=-1) #should use dark subtraction
        if not iterable(rep):
            print("reference failed")
            return
        fperc=percentile(rep,90,1)
        if any(fperc<10):
            message("No light for calibration!", title = 'Warning')
            return
        self.instr.flat=rep
        if rc.saturate>0:
            for i in range(len(self.instr.flat)):
                if any(self.instr.flat[i]>=rc.saturate):
                    self.display("channel %i saturated"%i)
        if self.combine and hasattr(self.instr,"setflat") and self.paren.spectrac.chansel<0:
            self.instr.setflat(perc=rc.flat_above_quantile,fresh=False)
            self.display("reference table calculated")
            if iterable(self.instr.ysel):
                self.display("channels :"+' '.join(["[%i]"%(sum(a) if iterable(a) else 0) for a in self.instr.ysel]))
                self.display("weights :"+' '.join(["[%.2f / %i]"%(a.sum(),sum(a.sum(0)>0)) for a in self.instr.transtable if iterable(a)]))
        else: 
            for i in range(len(self.instr.flat)):
                perc=percentile(self.instr.flat[i],[10,90])
                self.display("ch %i:%.0f-%.0f"%(i,perc[0],perc[1]))
            self.display("reference stored")
        self.paren.status_string="calibration finished"
        if not self.combine:
            for i in range(len(self.instr.flat)):
                if iterable(self.instr.ysel) and len(self.instr.ysel)>=i: 
                    self.instr.ysel[i]*=self.instr.flat[i]>=rc.underate
                self.instr.flat[i][self.instr.flat[i]<1]=1 #correct for zeros
        if self.config!=None: self.config.adjust_image()

intwid=50
class Scan(HasTraits):

    instr=None

    Xpts = Int(10, label="X points", desc="fast axis/radial")
    Ypts = Int(10, label="Y points", desc="slow axis/azimuth")
    Xstep = Float(rc.cart_step, label="X step", desc="instrumental steps [mm]")
    Ystep = Float(rc.cart_step, label="Y step", desc="instrumental steps [mm]")
    centerstart = Bool(True, label="start at center")
    callast = Bool(False, label="reference at last", desc="to measure reference when the scan is finished")
    zigzag = Bool(True, label="scan in both directions")
    radius = Float(50, label="wafer radius", desc="in [mm]")
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
    getpos = Button("Get Position")
    gopos = Button("Goto")

    up = Button("Up")
    down = Button("Down")
    left = Button("Left")
    right = Button("Right")
    near = Button("Near")
    far = Button("Far")
    shome = Button("Reset axis")
    gocenter = Button("Go center")
    srefer = Button("Calib sample")
    ptskip = Button("Skip")
    #res=[Button("%i s"%s) for s in [2,5,10]]
    view = View(VFold(
                HGroup(Item('Xpts',width=intwid),Item('Ypts',width=intwid),Item('Xstep',width=intwid),Item('Ystep',width=intwid)),
                #HGroup(Item('Xstep'),Item('Ystep')),
                HGroup(Item('radius', enabled_when="centerstart"),Item('centstr',style='readonly'),
                    Item('gocenter',show_label=False),Item('setcenter',show_label=False),enabled_when="ready"),
                HGroup(Item('centerstart'),Item('callast'),Item('zigzag'),Item('design',show_label=False),Item('dump',show_label=False),Item('pclear',show_label=False)),
                HGroup(Item('up',show_label=False),Item('down',show_label=False),Item('left',show_label=False),Item('right',show_label=False),
                    Item('near',show_label=False),Item('far',show_label=False),enabled_when="ready"),
                HGroup(Item('npoints',style='readonly'),Item('shome',show_label=False),Item('srefer',show_label=False),
                     Item('getpos',show_label=False),Item('cpos'),Item('gopos',show_label=False),enabled_when="ready"),
                    #Item('cmeasure',show_label=False)
                #Item('speed',show_label=True),
                label="Scanner", layout='split',
                )
            )
    program = []
    actpos = List(Float,[0,0])
    zpos = 0.
    count = [0,0]
    exelist = {}
    since_calib=0

    def setup(self):
        if self.instr==None:
            self.instr=labin.specscope()
            self.instr.adrive(rc.motor_lin_rate,rc.motor_lin_gap)# may not be present
        self.ready=True
        self.centpos=rc.xy_cent
        self.centstr=str(list(self.centpos))

    def _design_fired(self):
        from numpy import mgrid,newaxis,array
        #self.centpos=rc.xy_cent#self.setup()
        self.centstr=str(list(self.centpos))
        self.program=mgrid[:self.Xpts,:self.Ypts]*(array([self.Xstep,self.Ystep]).reshape(2,1,1))#.astype('float64')
        if self.zigzag:
            for i in range(1,self.Xpts,2):
                self.program[1][i]=self.program[1][i][::-1]
        center=array(self.centpos).reshape(2,1,1).astype('int32')
        if self.centerstart:
            self.program=self.program+center
            self.program[0]-=(self.Xpts*self.Xstep/2.)
            self.program[1]-=(self.Ypts*self.Ystep/2.)
        if self.radius>0:
            #for i in range(2):
            #    self.program[i]+=self.centpos[i]
            #if self.radius>0:
            #self.program=array(self.program)
            sel=((self.program-center)**2).sum(0)<self.radius**2
            print("selected %i points"%sum(sel))
            self.program=self.program[:,sel]
        else:
            self.program=self.program.reshape(2,self.Xpts*self.Ypts) #zrusena transpozice
        sel=self.program[0]>0
        sel*=self.program[1]>0
        sel*=self.program[0]<rc.xy_size[0]
        sel*=self.program[1]<rc.xy_size[1]
        self.program=list(self.program[:,sel].T)
        if self.callast:
            self.program.append([int(a) for a in rc.refer_pos[:2]])
        self.npoints=len(self.program)
        
        if 'plan' in self.exelist: #vykreslit
            from scanner import reband
            self.exelist['plan'].experiment.display("%i points out of accessible area"%(len(sel)-sum(sel)))
            print("plotting %i points (%i removed)"%(len(self.program),len(sel)-sum(sel)))
            self.exelist['plan'].design_show(self.program)
            self.exelist['plan'].experiment.clear_stack()
            self.exelist['plan'].wafer=reband.Wafer("",[])
        self.since_calib=0

    def disp_pos(self):
        #self.cpos=str([int(a) for a in list(self.actpos)])
        from numpy import round
        self.cpos=str([round(float(a),1) for a in list(self.actpos)])

    def _srefer_fired(self):
        #self.instr.goto(rc.refer_pos[0],rc.refer_pos[1]) #to be tested!!
        from numpy import array
        self.instr.goto(rc.refer_pos[0],rc.refer_pos[1],rc.refer_pos[2])#_gopos_fired()
        self.actpos=rc.refer_pos[:2]
        self.disp_pos()
        #self.instr.goto(rc.refer_pos[0],rc.refer_pos[1],self.instr.gz)

    def _gopos_fired(self):#_cpos_changed(self):
        #modified center position
        gpos=self.cpos.replace('[','').replace(']','').strip().split(',')
        if len(gpos)==2:
            try:
                pos=[int(float(g)) for g in gpos]
                self.instr.goto(pos[0],pos[1])
            except:
                pass

    def _shome_fired(self):
        # reset axis
        #for i in range(2):
        #    self.actpos[i]=rc.xy_cent[i]
        from numpy import iterable
        if self.instr==None: return 
        if iterable(rc.refer_pos):
            self.instr.ahome(rc.refer_pos)
            self.actpos=rc.refer_pos[:2]
        else:
            self.instr.ahome()
            self.actpos=rc.xy_cent
        self.disp_pos()
        self.setup()

    def _getpos_fired(self):
        self.instr.awrite("M114")
        inp=self.instr.acomm()
        for l in inp:
            sl=str(l)
            if sl.find('X')>0:
                #print("got resp "+sl)
                for k in sl.split('.0'):
                    if k.find('X')>=0:
                        ival=float(k[k.find(':')+1:])
                        self.actpos[0]=float(ival)
                    elif k.find('Y')>=0:
                        ival=float(k[k.find(':')+1:])
                        self.actpos[1]=float(ival)
                self.disp_pos()

    def _cmeasure_fired(self):
        print("trying to find wafer center [currently %s]"%str(self.centpos))
        self.centpos=self.actpos
        self.centstr=str(list(self.centpos))

    def _setcenter_fired(self):
        print("current position as wafer center")
        self._getpos_fired()
        if len(self.actpos)==2: #sanity check
            self.centpos=self.actpos.copy()
            self.centstr=str(list(self.centpos))
            print("center now at "+self.centstr)
        else:
            print("actual position misplaced - skipping")

    def _actpos_changed(self):
        self.disp_pos()

    def _gocenter_fired(self):
        from numpy import array
        self.instr.goto(self.centpos[0],self.centpos[1])
        self.actpos=self.centpos
        self.disp_pos()

    def _right_fired(self):
        if self.instr: self.instr.rate(1,self.Xstep)
        self.actpos[0]+=self.Xstep
        self.disp_pos()

    def _left_fired(self):
        if self.instr: self.instr.rate(1,-self.Xstep)
        self.actpos[0]-=self.Xstep
        self.disp_pos()

    def _down_fired(self):
        if self.instr: self.instr.rate(2,self.Ystep)
        self.actpos[1]+=self.Ystep
        self.disp_pos()

    def _up_fired(self):
        if self.instr: self.instr.rate(2,-self.Ystep)
        self.actpos[1]-=self.Ystep
        self.disp_pos()

    def _near_fired(self):
        if self.instr: 
            self.zpos-=rc.zaxis_step
            if self.zpos<0: self.zpos=0
            self.instr.awrite("G1 Z%.1f"%self.zpos)

    def _far_fired(self):
        if self.instr: 
            self.zpos+=rc.zaxis_step
            self.instr.awrite("G1 Z%.1f"%self.zpos)

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
        exper=self.exelist['plan'].experiment
        if len(self.program)==0: #end of plan
            if self.since_calib>0: #at least some points measured
                # will stop acquisition and save data
                self.exelist['plan']._start_stop_acquisition_fired()
                if exper.sname!="": exper._saveall_fired()
                #message(" scan finished! ", title = 'User mess')
                self.since_calib=0
            return
        if self.since_calib>0 and rc.save_period>0 and self.since_calib%rc.save_period==0:
            if exper.sname!="": exper._saveall_fired()
        point=self.program.pop(0)
        #movement
        if hasattr(self.instr,'ard'): #arduino/velleman
            if hasattr(self.instr,'goto'):
                self.instr.goto(point[0],point[1])
            else:
                self.instr.rate(1,point[0]-self.actpos[0])
                self.instr.rate(2,point[1]-self.actpos[1])
        else:
            print("riding %i,%i"%(point[0],point[1]))
            sleep(10-self.speed)
        #recalibration
        if exper.recalper>0 and self.since_calib>=exper.recalper:
            # recalibration planned
            if rc.refer_pos[0]*rc.refer_pos[1]>0: #have calibration sample
                self._srefer_fired() # goto reference sample
                sleep(5-self.speed)
                exper._refer_fired(interrupt=False) #measure reference
                self.since_calib=0        
        self.actpos=[float(p) for p in point]
        self.disp_pos()
        self.npoints-=1
        self.since_calib+=1
        #self.cpos=str(list(self.actpos))
        # run acquisition for all "keys" !!??
        for k in self.exelist.keys():
            self.exelist[k].acquire=True

    def height_calib(self,rep=1):
        rep=[]
        while self.npoints>0:
            self.next_point()
            rep.append(self.instr.focus(rep=rep,factval=[]))


#   def run_scan(self):
#        #obsolete
#        #now part of acquisition thread
#        self.process=self.instr.scan_save(self.Xstep,self.Ystep,self.Xpts,self.Ypts,self.radius)
#----------------------------------------

def analyse(func,odict,disp):
    """ Function called to do the processing (more time consuming, doesn't block the acquisition) - nothing implemented yet"""
    from scanner import spectra
    exec(func,odict)
    if 'samp' in odict: 
        disp("sample %s calculated"%odict['samp'].get_name())
    res=odict['results']
    if 'out' in res:
        output=str(res['out'])
        disp(output)
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
        self.analyse.calculate = True #call the Analyse event
        instr=self.experiment.instr
            
        if self.analyse.func!=None and self.analyse.bkgproc:
            odict=self.analyse.evalme(spect,prep_only=True)#{'results':{}}
            #print("starting proc. sample "+odict['samp'].get_name())
            #if instr!=None: odict['samp']=getattr(instr,"samp",None)
            self.processing_job = Thread(target=analyse, args=(self.analyse.func,odict,self.display))
            self.processing_job.start()
        else:
            self.processing_job = Thread(target=analyse, args=(spect,{},self.display))#,self.results))
            # do nothing for the moment
        self.result()

    def run(self):
        """ Runs the acquisition loop. """
        self.display('Spectrac started')
        self.n_img = 0
        #from numpy import array,iterable
        while not self.wants_abort:
            self.n_img += 1
            if hasattr(self,"prepare"): self.prepare() #next point
            spect=self.acquire(self.experiment)
            if not iterable(spect):
                print("acquisition failed")
                break #returned none
            #spect=array(spect)
            if self.experiment.refermat!='unity':
                instr=self.experiment.instr
                if instr.samp!=None:
                    instr.samp.calib()
                    if iterable(instr.transtable):
                        instr.samp.prof=instr.samp.trans(instr.transtable,instr.ysel,lev=rc.band_minlev,ref=instr.flat)
                    spect=instr.samp.prof.vals
                else:
                    print("sample not created")
            self.display('%d spectra captured' % self.n_img)
            if self.experiment.combine==False:
                if sum(array(self.experiment.instr.intfact)>0)==1: #just one channel
                    isel=[i for i in range(len(spect)) if self.experiment.instr.intfact[i]>0][0]
                    spect=spect[isel] 
            if hasattr(self,"stack"):
                sname=self.experiment.stack_name()
                if not sname in self.stack:
                    self.experiment.paren.stack_list.append(sname)
                self.stack[sname]=spect
            self.image_show(spect)
            self.process(spect)
        self.display('Spectrac stopped')

#--------------------------------------------------------------------
class ChanHandler(Handler):
    channels=List(Str)
    def object_instr_changed(self, info):
        self.channels=['all channels']+info.object.chanlist
        info.object.singlechan=self.channels[0]
    def object_ready_changed(self, info):
        self.channels=['all channels']+info.object.chanlist
        info.object.singlechan=self.channels[0]

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
    setmeup = Button("Initialize")
    setmeall = Button("Init and calibrate")
    equalize = Button("Equalize channels")
    calib = Button("Calibrate SiO2")
    focus = Button("Get focus")
    singlechan=Str#,Enum(['all channels','single channel'])
    chanlist=['single channel']
    chanshow = Button("Show channels")
    resol = Float(rc.meas_comb_resol, label="spectral resolution in eV")
    cooler=Bool(False)
    cooltemp = Range(-20,10,-10,label="Cooling",desc="required temp. in centigrade")
    norm = Tuple(Float,Float,Float)#List(Float, [1.,1,1])
    darkcorr = Bool(True, label="Dynamic dark")
    nonlincorr = Bool(True, label="Nonlin. correction")
    ready=Bool(False)
    chanmatch = Button("Match")
    chansel = -1
    refer_x = Int(rc.refer_pos[0])
    refer_y = Int(rc.refer_pos[1])
    cooler=Bool(False)
    debug=Bool(False)
    relative=Bool(rc.relative,label="Relative value")
    #norm = Array
    view = View(VGroup(
                Item('erange', editor=BoundsEditor(low_name = 'elow', high_name = 'ehigh')),
                Item('simulate'),
                HSplit(Item('resol',width=5),Item('cooltemp',width=5, enabled_when='cooler==True'),),
                HSplit(Item('darkcorr',width=5),Item('nonlincorr',width=5),),
                HSplit(Item('norm',style='simple',label="channel equalizer", enabled_when='instr!=None'),
                    #Item('singlechan',show_label=False, editor=ListEditor(style='readonly'),style='simple', enabled_when='instr!=None'),
                    Item('singlechan',show_label=False,editor=EnumEditor(name='handler.channels')),
                    Item('equalize',width=5,show_label=False),
                    Item('chanshow',width=5,show_label=False),Item('chanmatch')),#,enabled_when="instr.flat!=None")),
                Item('calib',width=5,show_label=False, enabled_when='instr.samp!=None'),Item('focus',show_label=False),
                HSplit(Item('refer_x'),Item('refer_y'),),
                HSplit(Item('debug'),Item('relative'),),
                HSplit(Item('setmeup',show_label=False),Item('setmeall',show_label=False)),#Item('port'),
                springy=True),
                resizable=True,
                handler=ChanHandler
            )#, enabled_when="config!=None")
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
                #self.instr=labin3d.webocean3d(path="http://localhost:%i/?exp="%self.port,chan=rc.chan_sel)
                self.instr=labin3d.ocean3d()
                if self.instr.dsize<=0:
                    message("Turn on the spectrometer, please ", title = 'User request')
                    self.instr=None 
                    return
                self.instr.gxmax,self.instr.gymax=rc.xy_size
            else: 
                #self.instr=labin.uniocean(path="http://localhost:%i/?exp="%self.port,chan=rc.chan_sel) ##FM replaced webocean
                self.instr=labin.oceanjaz()
        self.instr.setup([self.elow,self.ehigh],estep=self.resol,integ=self.exper.expo,aver=self.exper.aver)
        import os
        if hasattr(rc,'dark_path') and os.path.exists(rc.dark_path): #precalibration
            from numpy import loadtxt
            self.instr.dark=loadtxt(rc.dark_path,unpack=True)
            print('preloaded dark '+str(self.instr.dark.mean(1)))
        if len(self.instr.chanene)>0: self.exper.display(" found %i channels ..."%len(self.instr.chanene))
        #if len(self.instr.chanrange)>0: self.exper.display(" found %i channels ..."%len(self.instr.chanrange))
        self.exper.instr=self.instr
        if hasattr(self.instr,'ard'):
            self.exper.display("scanner connected to %s"%rc.ard_port)
            self.paren.scanner.instr=self.instr
            self.paren.scanner.ready=True
        self.pixels=self.instr.pixtable.copy()#
        #if self.instr.config.m_StartPixel>0 or (self.instr.config.m_StopPixel>1 and self.instr.config.m_StopPixel<len(self.pixels)):
        #    self.pixels=self.pixels[self.instr.config.m_StartPixel:self.instr.config.m_StopPixel]
        #self.instr.config.Material=b'simu'
        self.instr.samp=None
        self.exper.display("setup ok [%i pixels] ..."%len(self.pixels))
        self.norm=rc.chan_equal
        self.paren.status_string="now calibrate [Experiment/Dark + Reference]"
        self.chanlist=[('%.2f-%.2f eV'%(c.min(),c.max())) for c in self.instr.chanene]
        self.ready=True
        #for c in self.chanlist:
        #    self.singlechan.append(c)

    def _setmeup_fired(self):
        self.setup()
        self.exper.ready=True
        self.paren.ready=True
        #self.scanner.setup()

    def _setmeall_fired(self):
        self.setup()
        self.exper.ready=True
        self.paren.ready=True
        self.paren.scanner._shome_fired()
        self.exper._refer_fired()
        self.paren.scanner._gocenter_fired()
        #self.scanner.setup()

    def _focus_fired(self):
        rep=labin3d.zcalib(self.instr,rc.refer_pos[:2])
        print(rep[0])
        mid=(rep[0][0]+rep[0][1])/2.
        if mid<self.instr.gzmin: return
        if mid>self.instr.gzmax: return
        self.goto(rc.refer_pos[:2]+[mid])
        self.instr.gz=rc.refer_pos[-1]
        self.awrite("G92 Z%.2f"%rc.refer_pos[-1]) #set new position

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
        if self.instr!=None: 
            self.instr.intfact=list(self.norm)
            self.instr.set_acquisition(self.exper.expo,self.exper.aver)
 
    def _debug_changed(self):
        if self.instr!=None: 
            labin.loud=1 if self.debug else 0
 
    def _relative_changed(self):
        rc.relative=self.relative
 
    def _singlechan_changed(self):
        from numpy import ones
        if self.singlechan=='all channels': 
            self.norm=tuple(list(ones(len(self.norm))))
            self.chansel=-1
        else:
            isel=self.chanlist.index(self.singlechan)
            if isel<0: isel=len(self.norm)-1
            self.chansel=isel
            self.norm=tuple([1. if i==isel else 0 for i in range(len(self.norm))])
            if self.instr!=None: 
                self.instr.pixtable=self.instr.chanene[isel]
                self.pixels=self.instr.pixtable.copy()
                self.exper.display("chanel %i selected [%i pixels]"%(isel,len(self.pixels)))
        #self._norm_changed()

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
        if hasattr(self.instr,'transtable') and iterable(self.instr.transtable):
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
        thk=self.instr.samp.bands[0].guess()
        print("thickness:"+str(thk))
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
                allspect=[]
                for i in range(experiment.median):
                    allspect.append(self.instr.result(div_flat=rc.relative))
                spect=median(allspect,0)
                self.paren.display("Median")
            else:
                spect=self.instr.result(div_flat=rc.relative)
            if not(hasattr(spect,"shape")): return None
            if len(spect.shape)>1:
                return spect.sum(axis=0)
        else:
            if experiment.median>1:
                from numpy import median
                allspect=[]
                for i in range(experiment.median):
                    allspect.append(self.instr.result(div_flat=rc.relative))
                    experiment.display("Collecting %i/%i"%(i+1,experiment.median))
                spect=median(allspect,0)
                #self.paren.display("Median")
            else:
                spect=self.instr.result(div_flat=rc.relative,maxit=-1,smooth=experiment.smooth)
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
        self.instr.samp=self.instr.makesamp(spect)
        self.instr.samp.laystruct=experiment.material
        self.instr.samp.pos=list(self.paren.scanner.actpos.copy())
        if hasattr(self.paren,'wafer'): 
            self.paren.wafer.samps.append(self.instr.samp)
            self.instr.samp.wafer=self.paren.wafer
            self.instr.samp.update_nearest()
        return spect

    def get_match(self, caltab):
        perc=rc.flat_above_quantile
        ou,ol=labin.gettrans(self.instr.pixtable,self.instr.chanene,self.instr.flat,skiplowest=perc,rep=-1)
        data=self.instr.result(div_flat=rc.relative,maxit=-2) # reuse last measurement
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
    chuse = Enum([1,2,3], label="Channel to use")
    smnum = Int(0, label="Sample number to use")
    doplot = Bool(False, label="Show graph")
    bkgproc = Bool(False, label="Process", desc="in background")
    code = Code("return [xdat.mean(),ydat.mean()]")
    run = Button('Eval')
    debug = Button('Debug')
    res = Button('Reset')
    fitshow = Button('Show fit')
    variab = String('thick',label='What to show')
    fname = File('',label="Script name")
    load = Button('Load')
    save = Button('Save')
    dump = Button('Dump', desc="save displayed map to file")
    map = Button('Map')
    output = String("here comes the data",label="Output")
    view = View(
                HGroup(Item('erange', editor=BoundsEditor(low_name = 'elow', high_name = 'ehigh')),
                    Item('chuse'),Item('smnum'),Item('doplot'),Item('bkgproc')),
                Item('code',show_label=False),
                HGroup(Item('run',show_label=False),Item('res',show_label=False),Item('debug',show_label=False),
                    Item('map',show_label=False),Item('fitshow',show_label=False),Item('variab')),
                HGroup(Item('fname',show_label=False),Item('load',show_label=False),Item('save',show_label=False),Item('dump',show_label=False),),
                Item('output',show_label=False,style='readonly')
                )#, enabled_when="config!=None")
    calculate = Event
    func=None
    paren=None
    vals=[]

    def evalme(self,ydat=None,prep_only=False): 
        '''not working for multichannel
        '''
        samp=None
        if not iterable(ydat): ydat=self.instr.last
        if hasattr(self,'instr'):
            if len(self.instr.chanene)>0:
                if self.chuse>0:
                    if self.chuse>len(self.instr.chanene): self.chuse=1
                    xdat=self.instr.chanene[self.chuse-1]
                    ydat=ydat[self.chuse-1]
            else:
                xdat=self.instr.pixtable
        else: #processing loaded data
            if hasattr(self.paren,'wafer') and not iterable(ydat):
                if self.smnum>=len(self.paren.wafer.samps): self.smnum=0
                samp=self.paren.wafer.samps[self.smnum]
                if self.chuse>len(samp.bands): self.chuse=1
                if self.chuse>0: xdat=samp.bands[self.chuse-1].ix
                if iterable(ydat): 
                    if len(ydat.shape)==2: ydat=ydat[self.chuse-1]
                else: ydat=samp.bands[self.chuse-1].iy

        sel=(xdat>=self.elow)*(xdat<=self.ehigh)
        if not iterable(ydat) or len(ydat)!=len(sel): 
            if prep_only: return {}
            else: return 0
        odict={'xdat':xdat[sel],'ydat':ydat[sel],'samp':samp,'results':{}}
        if samp==None and hasattr(self.instr,"samp"): odict['samp']=self.instr.samp
        if prep_only: return odict
        exec(self.func,odict)
        return odict['results']

    def _reset_fired(self):
        self.vals=[]

    def _load_fired(self):
        import os
        if os.path.exists(self.fname):
            self.code=open(self.fname).read()
            self.status_string="Code loaded"
        else:
            self.status_string="File not found"
        
        
    def _save_fired(self):
        ofile=open(self.fname,"wb")
        ofile.write(self.code.encode("u8"))
        self.status_string="Code saved"
        ofile.close()

    def _debug_fired(self):
        try:
            self.output=str(eval(self.code))
        except:
            self.status_string="Error in expression"
        self.vals=[]

    def _fitshow_fired(self):
        if self.instr==None: return
        if hasattr(self.instr,"samp"): 
            samp=self.instr.samp
        else:
            if not(hasattr(self.paren,'wafer')) or self.paren.wafer==None: return 
            if self.smnum>=len(self.paren.wafer.samps): self.smnum=0
            samp=self.paren.wafer.samps[self.smnum]
            
        #self.paren.figure.clean()
        ax=self.paren.figure.axes[0]
        ax.clear()
        #self.paren.spect_last=[]
        if len(samp.thick)==0: 
            # no fit yet
            self.status_string="No fit saved yet!"
            return
        rep=samp.plot(ax=ax)
        self.paren.spect_last=rep
        self.paren.figure.canvas.draw()


    def _run_fired(self):
        #rep=exec(self.code) not working
        imp=self.code
        imp=imp.replace("return ","results['out']=")
        self.func=compile(imp,'compiled','exec')
        self.output=str(self.evalme()['out'])
        #self.output="Evaluation failed"

    def _map_fired(self):
        '''apply calculation on all data
            variab: precalculated values to display (from panel.wafer)
            code: apply code to all 
        '''
        self.vpos=[]
        self.vals=[]
        import numpy as np
        if len(self.variab)>0 and hasattr(self.paren,'wafer'): #already calculated values
            varfunc=compile("result="+self.variab,'compiled','exec')
            print("calculating '%s' for %i points"%(self.variab,len(self.paren.wafer.samps)))
            for samp in self.paren.wafer.samps:
                if self.variab=='thick':
                    if len(samp.thick)<1: continue
                    zval=np.mean(list(samp.thick.values()))
                else:
                    odict={'result':0,'bands':samp.bands}
                    exec(varfunc,odict)
                    zval=odict['result']
                    if np.iterable(zval): zval=zval[0]
                if len(samp.pos)<2: continue
                #print("%i values calculated"%(len(self.vals)))
                self.vpos.append(samp.pos)
                self.vals.append(zval)
        else:
            imp=self.code
            imp=imp.replace("return ","results['out']=")
            self.func=compile(imp,'compiled','exec')
            odict={'results':{}}
            if hasattr(self.paren,'wafer'): #using wafer 
                for samp in self.paren.wafer.samps:
                    odict['samp']=samp
                    exec(self.func,odict)
                    self.vpos.append(samp.pos)
                    self.vals.append(odict['out'])
            elif len(self.paren.experiment.stack)>0: #using stack
                for n in self.paren.experiment.stack.keys():
                    spec=self.paren.experiment.stack[n]
                    try:
                        pos=[int(float(i)) for i in n.split('_')[1:]]
                        self.vpos.append(pos[:2])
                    except:
                        continue
                    odict['ydat']=spec
                    exec(self.func,odict)
                    self.vals.append(odict['out'])
        
        print("showing values from %.1f to %.1f"%(np.min(self.vals),np.max(self.vals)))
        self.paren.design_show(self.vpos,self.vals)
    
    def _dump_fired(self):
        from numpy import array,savetxt,concatenate
        if len(self.vpos)==len(self.vals):
            odata=concatenate([array(self.vpos).T,[self.vals]])
        else:
            odata=self.vals
        savetxt(self.fname,odata.T,fmt="%8.3f")

    def _calculate_fired(self):
        '''Event called by analyse
        odict: input data
        is processed in run thread - not for time consuming analysis...
        '''
        #if self.func!=None:
        if False:
            odict={'results':{'out':[]}}
            odict['samp']=getattr(self.instr,"samp",None)
            exec(self.func,odict)
            res=odict['results']
            if 'out' in res:
                self.output=str(res['out'])

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
    status_string = String()
    result_plot=[]
    grcolors="rgbkmc"
    showme=Button("Show")
    delme=Button("Delete")
    loadme=Button("Load set")
    sname = File(rc.test_map,label="Filename")
    saveme=Button("Save selected")
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
                        HGroup(Item("showme", show_label=False),Item("delme", show_label=False)),
                        HGroup(Item("saveme", show_label=False),Item("loadme", show_label=False),
                            Item("sname", show_label=False,editor = FileEditor(dialog_style='load'))),
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
        self.analyse.paren = self
        # self.scanner.exelist['spectrac']=self.spectrac 
        # only for automatic acquisition
        self.scanner.exelist['plan']=self
        #self.bardial.open()

    def _showme_fired(self):
        '''
        show spectra from stack
        '''
        from matplotlib import patches
        for p in self.stack_disp:
            try:
                p.remove()
                del p
            except:
                pass
        self.stack_disp=[]
        i=0
        dax=self.design.axes[0]
        for a in self.stack_selec:
            if not(a) in self.experiment.stack: continue
            print(a,sum(self.experiment.stack[a][0]))
            if self.experiment.combine:
                mcol=rc.all_color[i%len(rc.all_color)]
                p=self.figure.axes[0].plot(self.spectrac.pixels,self.experiment.stack[a],mcol)[0]
                self.stack_disp.append(p)
            else:
                mcol=rc.all_color[i%len(rc.all_color)]
                for j in range(len(self.experiment.instr.chanene)):
                    ix=self.experiment.instr.chanene[j]
                    if self.experiment.instr.intfact[j]==0: continue
                    p=self.figure.axes[0].plot(ix,self.experiment.stack[a][j],mcol)[0]
                    self.stack_disp.append(p)
            coors=a.split("_")
            if len(coors)==3:
                try:
                    qx,qy=float(coors[1]),float(coors[2])
                    p=patches.CirclePolygon([qx,qy],1,color=mcol,alpha=0.3)
                    dax.add_patch(p)
                    self.stack_disp.append(p)
                except:
                    print("failed parse of "+a)
                    pass
            i+=1
        
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

    def _loadme_fired(self):
        from scanner import reband
        import os
        if not os.path.exists(self.sname):
           message("file "+self.sname+" not found")
           return
        qname=str(self.sname)
        import re
        ext=qname[qname.rfind(os.path.extsep):]
        qname=re.sub("_\d"+ext,"_%i"+ext,qname)
        elist=[]
        for i in range(5):
            if os.path.exists(qname%i):
                elist.append(i)
        print("loading from "+qname+str(elist))
        self.wafer=reband.Wafer(qname,elist,laystruct=self.experiment.material,maxband=-1,headerow=2,position=2)
        sm=self.wafer.samps[len(self.wafer.samps)//2]
        sm.calib()
        if self.experiment.instr==None:
            self.experiment.instr=labin.specscope()
            self.experiment.instr.chanene=[bd.ix for bd in sm.bands]
            self.experiment.instr.intfact=list(self.spectrac.norm)
        self.wafer.make_hash()
        bad=[]
        plan=[]

        ax=self.design.axes[0]
        def accept(event): #for lasso selection
            if event.key == "enter":
                print("Selected points:")
                print(self.selector.xys[self.selector.ind])
                self.selector.disconnect()
                ax.set_title("")
                ax.figure.canvas.draw()
            self.status_string="%i points selected"%(len(self.selector.ind))

        def onclick(event):
            from matplotlib import patches
            import numpy as np
            #global aw
            #self=aw.panel
            if not hasattr(self,'wafer'): return
            ssel=self.wafer.get_nearest(event.xdata, event.ydata)
            if len(ssel)==0: return
            samp=ssel[0]
            print('%s click: %s button=%d, xdata=%f, ydata=%f' %
                ('double' if event.dblclick else 'single', samp.get_name(), event.button,
                samp.pos[0], samp.pos[1]))
            dax=self.design.axes[0]
            ip=len(dax.patches)
            mcol=rc.all_color[ip%len(rc.all_color)]
            if event.button==1:
                if samp.get_name() in self.stack_selec:
                    self.stack_selec.remove(samp.get_name())
                    for k in range(len(self.stack_disp)):
                        p=self.stack_disp[k]
                        if type(p)==patches.CirclePolygon:
                            if np.all(p.xy==list(samp.pos)):
                                p.remove()
                                del self.stack_disp[k]
                                break
                else:
                    p=patches.CirclePolygon(samp.pos,1,color=mcol,alpha=0.3)
                    dax.add_patch(p)
                    self.stack_disp.append(p)
                    self.stack_selec.append(samp.get_name())
                ax.figure.canvas.draw()

        cid = self.design.canvas.mpl_connect('button_press_event', onclick)

        for sm in self.wafer.samps: 
            k=sm.get_name()
            data=[b.absol() for b in sm.bands]
            for i in range(len(data)):
                data[i][data[i]<0]=0 #rc.min_spec
            self.experiment.stack[k]=data
            self.stack_list.append(k)
            coors=k.split("_")
            if len(coors)==3:
                try:
                    qx,qy=float(coors[1]),float(coors[2])
                    plan.append([qx,qy])
                except:
                    bad.append(k)
        if len(plan)>0:
            #pts=ax.plot([p[0] for p in plan],[p[1] for p in plan],'bd-')[0]
            pts=ax.scatter([p[0] for p in plan],[p[1] for p in plan],marker='d')
            ax.figure.canvas.draw()
            #self.selector = SelectFromCollection(ax, pts)
            #self.design.canvas.mpl_connect("key_press_event", accept)

        self.experiment.display("%i points loaded"%(len(plan)))

    def _saveme_fired(self):
        self.experiment._saveall_fired(selected=True)

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
            if self.spect_last!=None: #cleaning
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

    def adjust_image(self,margin=0.02,percent=99):
        from numpy import iterable,percentile
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
                if percent>0 and percent<100:
                    ylow=percentile(yval,100-percent)
                    yhig=percentile(yval,percent)
                else:
                    ylow,yhig=yval.min(),yval.max()
                if ymin>ylow: ymin=ylow
                if ymax<yhig: ymax=yhig
            if margin>0:
                ydis=ymax-ymin
                ymin,ymax=ymin-ydis*margin,ymax+ydis*margin
            if ymin<rc.graph_min: ymin=rc.graph_min
            if rc.relative and (ymax>rc.graph_max): ymax=rc.graph_max
            self.figure.axes[0].set_ylim(ymin,ymax)
            self.figure.canvas.draw()

    def add_line(self, string):
        """ Adds a line to the textbox display."""
        self.results_string = (string + "\n" + self.results_string)[0:1000]

    def image_show(self, spect):
        """ Plots an image on the canvas in a thread safe way. """
        #self.figure.axes[0].images=[]
        if self.spect_last==None:
            pin=self.spectrac.instr
            if rc.debug>0: print("(re)new graphs")
            if type(spect)!=list and len(spect.shape)==1:#self.experiment.combine:
                if len(self.spectrac.pixels)==len(spect): pixs=self.spectrac.pixels
                else: 
                    ipck=arange(len(pin.chanene))[array(pin.intfact)>0]
                    pixs=pin.chanene[ipck[0]]
                self.spect_last=self.figure.axes[0].plot(pixs,spect,rc.line_color,lw=rc.line_width)[0]
            else:
                from numpy import array,arange
                self.spect_last=[]
                ipck=arange(len(pin.chanene))[array(pin.intfact)>0]
                ibnd=min(len(ipck),len(spect))
                for i in range(ibnd):
                    dx,dy=pin.chanene[ipck[i]],spect[i]
                    if pin.ysel!=None:
                        dx,dy=dx[pin.ysel[i]],dy[pin.ysel[i]]
                    self.spect_last.append(self.figure.axes[0].plot(dx,dy,lw=rc.line_width,color=self.grcolors[i])[0])
            self.figure.axes[0].grid(1)
        else:
            if type(spect)!=list and len(spect.shape)==1:#self.experiment.combine:
                self.spect_last.set_ydata(spect)
            else:
                from numpy import array
                pin=self.spectrac.instr
                nbnd=sum(array(pin.intfact)>0)
                if rc.debug>1: print("updating %i graphs"%nbnd)
                for i in range(len(spect)):
                    #if not iterable(pin.ysel[i]): continue  
                    if len(self.spect_last)<=i: break
                    if self.spectrac.norm[i]>0: #instr.intfact
                        if pin.ysel!=None:
                            self.spect_last[i].set_ydata(spect[i][pin.ysel[i]])
                        else:
                            self.spect_last[i].set_ydata(spect[i])
        self.figure.canvas.draw()
        #self.profile.canvas.draw()

    def result_update(self):
        """plots results of the analysis"""
        from numpy import arange,iterable
        if not(self.analyse.doplot): return
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

    def design_show(self, plan, values=[]):
        ax=self.design.axes[0]
        ax.clear()
        if len(values)>0:
            ax.figure.clear()
            sc=ax.scatter([p[0] for p in plan],[p[1] for p in plan],c=values,cmap=glob_cmap)
            print("scaling to "+str(self.scanner.xdim)+" and "+str(self.scanner.ydim))
            ax.figure.colorbar(sc)
        elif len(plan)>0:
            ax.plot([p[0] for p in plan],[p[1] for p in plan],'bd-')
            pp=ax.get_xlim()
            mdis=(pp[1]-pp[0])*0.02
            self.scanner.xdim=[pp[0]-mdis,pp[1]+mdis]
            pp=ax.get_ylim()
            if pp[0]<pp[1]: pp=pp[::-1]
            self.scanner.ydim=[pp[0]-mdis,pp[1]+mdis]
        if hasattr(self.scanner,'xdim'): ax.set_xlim(*self.scanner.xdim)
        if hasattr(self.scanner,'ydim'): ax.set_ylim(*self.scanner.ydim)
        self.design.canvas.draw()

    def _acquire_fired(self):
        # only shows the point on display
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


import sys
if __name__ == '__main__':
    aw=MainWindow()
    #if "--init" in sys.argv: rc.auto_init=True
    for v in sys.argv:
        if v.find("chan")>0:
            rc.chan_sel=int(v[v.find('=')+1:])
    aw.configure_traits(view=tabbedview)
    message(" press Setup first ", title = 'User request')
        

def onclick(event):
    global aw
    self=aw.panel
    if not hasattr(self,'wafer'): return
    ssel=waf.get_nearest(event.xdata, event.ydata)
    if len(ssel)==0: return
    samp=ssel[0]
    print('%s click: %s button=%d, xdata=%f, ydata=%f' %
        ('double' if event.dblclick else 'single', samp.get_name(), event.button,
        samp.pos[0], samp.pos[1]))
    dax=self.design.axes[0]
    ip=len(dax.patches)
    mcol=rc.all_color[ip%len(rc.all_color)]
    if event.button==1:
        p=patches.CirclePolygon(samp.pos,1,color=mcol,alpha=0.3)
        dax.add_patch(p)
        self.stack_disp.append(p)
    self.stack_selec.add(samp.get_name())
