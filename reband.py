import numpy as np

refer=None
diel={}
minmodval=0.05
refthk=0.5#23.7 #oxide thickness on reference sample in [nm]
rescale=1
#indir="C:/Users/Admin/Documents/Lab/MOCVD/lpcvd-calib/"
indir="C:/Users/Optik/Documents/Data/Calib/"
sinname="sin_nk_test.mat"#"sinx_cl2.mat"

def transpose(self,nbext=0,pos=[]):
    if nbext==0: nbext=len(self.samps)
    gx=[sm.bands[0].ix for sm in self.samps]
    sm=self.samps[0]
    nsamps=[]
    for j in range(1,len(sm.bands)):
        gy=[sm.bands[j].iy for sm in self.samps]
        sdata=[]
        for i in range(nbext):
            sdata.extend([gx[i],gy[i]])
        nsamps.append(Sample(None,data=sdata,maxband=nbext))
        nsamps.wafer=self
        if j<len(pos):
            nsamps[-1].pos=pos[j]
    return nsamps

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
    #for i in range(len(iend)):
    #    print(iend[i]-istart[i])
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

class Band():
    ix=None
    iy=None
    samp=None
    qfit=0

    def __init__(self,samp,loc=0,bord=10):
        self.samp=samp
        if loc<0:
            self.ix,self.iy=samp.data[0],samp.data[-loc]
        else:
            self.ix,self.iy=samp.data[loc],samp.data[loc+1]
        #self.sel=slice(10,-10)
        self.sel=self.ix>0
        self.sel[:bord]=False
        self.sel[-bord:]=False
        self.correct=""
        self.scale=1
        self.rat=[1,0]
        
    def len(self):
        return self.ix.size()

    def range(self,imin=0,imax=None):
        self.sel=(self.ix>=imin)
        if imax!=None:self.sel*=(self.ix<=imax)

    def bord(self,left=0,right=0):
        good=np.arange(len(self.ix))[self.sel]
        self.sel[good[:left]]=False
        if right>0: self.sel[good[-right:]]=False

    def corrdark(self,dark):
        #dark measured relatively
        if self.correct.find('d')>=0: return
        if np.any(dark>1):
            bad=np.where(dark>=1)[0]
            imid=len(self.iy)/2
            if sum(bad<imid)>0:
                imin=bad[bad<imid][-1]+1
                self.sel[:imin]=False
            if sum(bad>imid)>0:
                imax=bad[bad>imid][0]-1
                self.sel[imax:]=False
        self.iy[self.sel]=(self.iy-dark)[self.sel]/(1-dark[self.sel])
        self.correct+='d'

    def absol(self):
        return self.iy*self.samp.norm(self.ix)

    def guess(self,dlist=None,rep=1):
        #dlist in nanometers
        from scanner import profit
        epsi=[diel['ksio2'](self.ix),diel['ksi'](self.ix)]
        if not(np.iterable(dlist)): dlist=np.arange(20,300,5)
        nor=lambda d:((self.absol())[self.sel]/self.model([d],renow=True)[self.sel]).std()
        slist=[nor(d) for d in dlist]
        if rep==0: return slist
        return dlist[np.argmin(slist)]

    def renorm(self,minval=0.05,thick=None,ngrp=15.,loud=0):
        if thick==None:
            if len(self.samp.thick)>0:
                thick=self.samp.get_thick()#np.median(list(self.samp.thick.values()))
            else: return
        modval=self.model([thick,1,0],renow=False)
        osel=self.sel*(modval>minval)
        yval=self.absol()
        osel*=yval>minval
        if sum(osel)<3:
            print("too few points")
            osel=self.sel.copy()

        #robust init (outlier lookup)
        from scipy import ndimage as nd
        mmin=modval[osel].min()
        mstep=(modval[osel].max()-mmin)/ngrp
        modlab=((modval[osel]-mmin)/mstep).astype(int) #group labels
        meds=nd.median(yval[osel],modlab,range(max(modlab)-1))
        if loud>1:
            print(meds)
            print(modlab.min(),modlab.max(),sum(modlab==modlab.max()))
            print(nd.standard_deviation(yval[osel],modlab,range(max(modlab)-1)))
        vpos=mmin+mstep*np.r_[0.5:ngrp]
        vpos=vpos[:len(meds)]
        res0=np.polyfit(vpos,meds,1)
        dif=abs(modval[osel]*res0[0]+res0[1]-yval[osel])
        lim=np.percentile(dif,90)
        isel=np.r_[:len(osel)][osel][dif>lim] #indices where residual is too big
        osel[isel]=False

        res=np.polyfit(modval[osel],yval[osel],1)
        if res[0]<0:
            #print("neg.scale - rejected")
            res=[1,0]
        self.rat=res

    def trans(self,ttable,tsel,lev=0.2,qspline=10.,pixval=None):
        '''converts pixels to new bins defined by transform. ttable (resolution)
        calculates uncertainties from dispersions within bins

        if orig. points have uncertainties, trying to propagate them
        '''
        #self.wei=ttable.sum(0)
        #self.wei[self.sel==False]=0
        self.wei=(ttable*self.sel[tsel][:,np.newaxis]).sum(0)
        prof=self.absol().copy()
        prof[self.sel==False]=0
        zoo=prof[tsel].dot(ttable)
        zoo2=(prof[tsel]**2).dot(ttable)
        qsel=self.wei>lev
        self.wei[qsel==False]=0
        qmid=zoo[qsel]/self.wei[qsel]
        self.sig=np.sqrt(zoo2[qsel]/self.wei[qsel]-qmid**2)
        if np.iterable(pixval): bins=pixval[qsel]
        else: bins=np.r_[:sum(qsel)]
        if qspline>0:
            from scipy import interpolate as ip
            smoothed=ip.UnivariateSpline(bins,self.sig,s=qspline)
            self.sig=np.abs(smoothed(bins))
        from uncme import uarray
        self.mid=uarray(qmid,self.sig)
        return self.mid

    def trans2(self,ttable,ysel,ref,lev=0.1,hasnan=False):
        zlst=[]
        from uncme import uarray
        #prof=self.absol()#.copy()
        uref=uarray(ref,np.sqrt(ref))
        uprof=uarray(self.absol()*ref,abs(self.iy)*np.sqrt(ref))
        prof=(uprof*self.samp.norm(self.ix))/uref
        self.wei=(ttable*self.sel[ysel][:,np.newaxis]).sum(0)
        qsel=self.wei>lev
        #print(sum(qsel))
        self.wei[qsel==False]=0
        if hasnan: nsel=np.isnan(prof.errs()[ysel])==False
        rlen=len(ttable[0])
        for j in np.arange(rlen)[qsel]:#[ttable.sum(0)>0]:
            xsel=(ttable[:,j]>0)
            #if hasnan: qsel*=nsel
            ok=prof.nums[ysel][xsel]
            zlst.append(sum(ok*ttable[:,j][xsel])/self.wei[j])
        self.mid=uarray(zlst)

    def model(self,pars,renow=False):
        from . import profit
        if len(pars)<3:
            if renow: self.renorm(minmodval,pars[0])
            if hasattr(self,"rat"): p=[pars[0],self.rat[0],self.rat[1]]
            else: p=[pars,1,0]
        else:p=pars[:3]
        if self.samp.lay!=None:
            if self.samp.lay.find('SiN')>=0:
                if not 'ksin' in diel:
                    from scipy import interpolate as ip
                    tsin=np.loadtxt(indir+sinname,unpack=True,skiprows=3)
                    diel['ksin']=ip.interp1d(tsin[0],(tsin[1]+1j*tsin[2])**2)
                epsi=[diel['ksin'](self.ix),diel['ksi'](self.ix)]
            elif self.samp.lay.find('cau')>=0:
                if len(pars)>3:
                    if not 'cau' in diel: diel['cau']=lambda x:np.polyval(pars[3:][::-1],x**2)**2
                    epsi=[np.polyval(pars[3:][::-1],self.ix**2)**2,diel['ksi'](self.ix)]
                elif 'cau' in diel:
                    epsi=[diel['cau'](self.ix),diel['ksi'](self.ix)]
            else:
                epsi=[diel['ksio2'](self.ix),diel['ksi'](self.ix)]
        else:
            epsi=[diel['ksio2'](self.ix),diel['ksi'](self.ix)]
        return profit.plate(self.ix,epsi,[p[0]])*p[1]+p[2]

    def plot(self,amodel=False,match=True):
        from matplotlib import pyplot as pl
        if match:
            pl.plot(self.ix[self.sel],(self.absol()[self.sel]-self.rat[1])/self.rat[0])
        else:    
            pl.plot(self.ix[self.sel],self.absol()[self.sel])
        if amodel and len(self.samp.thick)>0:
            thick=self.samp.get_thick()
            if match:
                pl.plot(self.ix[self.sel],(self.model([thick])[self.sel]-self.rat[1])/self.rat[0])
            else:
                pl.plot(self.ix[self.sel],self.model([thick])[self.sel])

    def fit(self,inval=None,irat=[1.3,-0.15],save=None,refer=None):
        from scipy import optimize as op
        if refer!=None:
            resid=lambda p:sum((self.absol()[self.sel]-self.model(p)[self.sel])**2/refer.iy[self.sel]**2)
        else:
            resid=lambda p:sum((self.absol()[self.sel]-self.model(p)[self.sel])**2)
        if inval==None: return resid
        from numpy import iterable
        if iterable(inval):
            zpar=op.fmin(resid,[inval[0],irat[0],irat[1]]+list(inval[1:]),full_output=False,disp=False)
        else:
            zpar=op.fmin(resid,[inval,irat[0],irat[1]],full_output=False,disp=False)
        self.qfit=resid(zpar)
        if save!=None:
            self.samp.thick[save]=zpar[0]
            self.samp.chi2[save]=self.qfit
            self.rat=list(zpar[1:])
        return zpar

def overlaps(caltab):
    from numpy import sum,where
    cn=(sum([c.sum(0)>0 for c in caltab],0)-0.9).astype(int)
    istart=where(cn[1:]>cn[:-1])[0]
    iend=where(cn[1:]<cn[:-1])[0]
    return list(istart),list(iend)

class Sample():
    #bands=[] #to je spatne - sdilena
    data=None
    fname=""
    lay=None
    wafer=None
    pos=[]
    inval_thick=set()

    def __init__(self,fname,laystruct="SiO2/Si",delim=None,maxband=0,data=None,headerow=0,bord=10):
        self.bands=[]
        self.thick={}
        self.chi2={}
        self.fname=fname
        self.nearest={}
        if fname==None:
            assert np.iterable(data)
            self.data=data
        else:
            import os
            if not(os.path.exists(fname)):
                print("file "+fname+" not found")
                return
            self.data=np.loadtxt(self.fname,unpack=True,delimiter=delim,skiprows=headerow)
            if self.data.shape[0]>self.data.shape[1]: self.data=self.data.T
        if maxband==0: maxband=len(self.data)//2
        if maxband<0:
            for i in range(-1,maxband,-1):
                self.bands.append(Band(self,i,bord=bord))
        else:
            for i in range(0,maxband*2,2):
                self.bands.append(Band(self,i,bord=bord))
        self.norm=lambda x:1
        if laystruct!=None: self.lay=laystruct

    def calib(self):
        from scanner import spectra,profit
        import pickle,os
        #epssi=spectra.dbload("cSi_asp")
        #pickle.dump(epssi,open(indir+"si_eps_full.mat","w"))

        from scipy import interpolate as ip
        if not 'ksi' in diel:
            if not os.path.exists(indir+"si_eps_fulld.mat"):
                epssi=spectra.dbload("cSi_asp")
            else:
                epssi=pickle.load(open(indir+"si_eps_fulld.mat","rb"))
            diel['ksi']=ip.interp1d(epssi[0],epssi[1])
        if not 'ksio2' in diel:
            if not os.path.exists(indir+"sio2_palik_g.mat"):
                x=np.r_[0.5:6.5:0.01]#epssi[0]
                tsio2=[x,np.polyval(spectra.cau_sio2,x)]
            else:
                tsio2=np.loadtxt(indir+"sio2_palik_g.mat",unpack=True,skiprows=3)
            diel['ksio2']=ip.interp1d(tsio2[0],tsio2[1]**2)
        if refthk<1.: self.norm= lambda px: profit.reflect(diel['ksi'](px))
        else: self.norm= lambda px: profit.plate(px,[diel['ksio2'](px),diel['ksi'](px)],[refthk])

    def corrdark(self,dark):
        for i in range(len(self.bands)):
            if i>=len(dark.bands): break
            self.bands[i].corrdark(dark.bands[i].iy)
    
    def get_thick(self,select=None):
        if len(self.thick)==0: return None
        vlist=[]
        if select!=None:
            vlist=[v for k,v in self.thick.items() if k.find(select)>=0]
        elif len(self.inval_thick)>0:
            vlist=[v for k,v in self.thick.items() if k not in self.inval_thick]
        if len(vlist)==0:
            vlist=list(self.thick.values())
        vlist=[v for v in vlist if v>0]
        return np.median(vlist)
            
    def get_qfit(self,thk=None,save=True):
        res=[]
        if thk==None: 
            thk=self.get_thick()
            if thk==None: return
        for bd in self.bands:
            if not hasattr(bd,'rat'):
                bd.renorm(thick=thk)
            bresid=bd.fit(None)
            res.append(bresid([thk]+list(bd.rat)))
            if save: bd.qfit=res[-1]
        return res
    
    def renorm(self):
        for bd in self.bands:
            bd.renorm()
    
    def fit(self,inval=None,save=None,refer=None):
        from scipy import optimize as op
        if refer!=None:
            bresid=[self.bands[i].fit(None,refer=refer.bands[i]) for i in range(len(self.bands))]
        else:
            bresid=[self.bands[i].fit(None) for i in range(len(self.bands))]
        resid=lambda p:sum([bresid[i]([p[0],p[2*i+1],p[2*i+2]]) for i in range(len(self.bands))])
        inarr=np.concatenate([[inval]]+[b.rat for b in self.bands])
        if inval==None: return inarr
        zpar=op.fmin(resid,inarr,full_output=False,disp=False)
        if np.isnan(zpar[0]):
            return zpar
        if save!=None:
                self.thick[save]=zpar[0]
        for i,b in enumerate(self.bands):
            p=zpar
            b.qfit=bresid[i]([p[0],p[2*i+1],p[2*i+2]])
            if save!=None:
                b.rat=[p[2*i+1],p[2*i+2]]
        return zpar

    def trans(self,ttable,tsel,lev=0.2,qspline=10.,pixval=None,reweig=True,ref=None):
        if np.iterable(ref):
            data=[b.trans2(ttable[i],tsel[i],lev=lev,ref=ref[i]) for i,b in enumerate(self.bands)]
        else:
            data=[b.trans(ttable[i],tsel[i],lev=lev,qspline=qspline,pixval=pixval) for i,b in enumerate(self.bands)]
        zelen=ttable[0].shape[1]
        from uncme import uarray
        res=uarray(np.zeros(zelen),np.zeros(zelen))
        ob=np.sum([b.wei for b in self.bands],0)
        #print(sum(ob<1),sum(ob==1))
        ob[ob<=0]=1
        for b in self.bands:#range(len(ttable)):
            b.wei/=ob
            nadd=b.mid.nums.copy()*b.scale
            if reweig: nadd*=b.wei[b.wei>lev]
            res.nums[b.wei>lev]=res.nums[b.wei>lev]+nadd
        return res

    def mismat(self,ttable,ret=0):
        istart,iend=overlaps(ttable)
        from uncme import uarray
        rep=[]
        #zelen=len(pixval)
        ival=np.r_[:len(ttable[0][0])]
        for i in range(len(istart)):
            #res=uarray(np.zeros(zelen),np.zeros(zelen))
            #xint=pixval[istart[i]]
            msel=(np.sum([(b.wei>0) for b in self.bands],0)>1)
            sel=(ival>istart[i])*(ival<iend[i])
            laps=[]
            for b in self.bands:
                #sel=(pixval>xint[0])*(pixval<xint[1])
                tsel=(b.wei>0.)
                if sum(tsel*sel)==0: continue
                laps.append(b.mid.nums[sel[tsel]*msel[tsel]])
                #rep.append([b.mid[istart[i]:iend[i]] for b in self.bands])
            if len(laps)<2: continue #should not happen
            if ret==1: rep.append((laps[1]/laps[0]))
            else: rep.append(uarray(laps[1]/laps[0]).mean())
        return rep

    def save(self,fname,liner=False):
        if liner:
            of=open(fname,"w")
            for i,b in enumerate(self.bands):
                of.write("# band %i"%i)
                for j in range(len(b.ix)):
                    of.write("%.3f  %.3f \n"%(b.ix[i],b.iy[i]))
            of.close()
        else:
            blen=max([len(b.ix) for b in self.bands])
            #self.data=np.concatenate([[self.bands.ix[i],self.bands.iy[i]] for i in range(len(self.chanene))]
            np.savetxt(self.data)

    def get_name(self):
        return "pos_"+str(list(self.pos))[1:-1].replace(",","_").replace(" ","")
            

    def census(self,lab='first',outliers=20,valid=[]):
        from numpy import percentile
        rvals=[sq.thick[lab] for sq in self.nearest.values() if hasattr(sq,"thick") and lab in sq.thick and sq.thick[lab]>0]
        if len(valid)>1:
            rvals=[t for t in rvals if t>=valid[0] and t<=valid[1]]
        if len(rvals)>2 and outliers>0:
            zmin,zmax=np.percentile(rvals,outliers),np.percentile(rvals,100-outliers)
            zmin,zmax=zmin*1.5-zmax*0.5,zmax*1.5-zmin*0.5
            for i in range(len(rvals)-1,-1,-1):
                if rvals[i]>zmax or rvals[i]<zmin:
                    del rvals[i]
        return rvals

        
    def dump_hash(self,ofile,label=None):
        if label==None and len(self.pos)>0:
            label=self.get_name()
        for k in self.thick.keys():
            if k in self.inval_thick: continue
            ofile.write(label+": %s : %.4f\n"%(k,self.thick[k]))
            
    def load(self,fname):
        of=open(fname)
        for l in of:
            if l[0]=="#" and l.find('band')>0:
                bid=int(l[l.rfind(' '):])
                while bid>=len(self.bands):
                    self.bands.append(Band())
            else:
                dat=[float(q) for q in l.split()]
        of.close()

    def plot(self,amodel=False,lims=[0,1],unit='eV',match=True):
        for b in self.bands:
            b.plot(amodel,match=match)
        from matplotlib import pyplot as pl
        pl.ylim(*lims)
        pl.xlabel(unit)
        pl.grid()

posfun=lambda ix,iy:np.any([(edgpos==[ix,iy]).sum(1)==2])

class Wafer():

    expos=0

    def __init__(self,pattern,irange,laystruct=None,delim=None,maxband=0,headerow=0,position=0):
        self.names={}
        import os
        self.samps=[Sample(pattern%i,laystruct=laystruct,delim=delim,maxband=maxband,headerow=headerow) for i in irange]
        self.samps=[sm for sm in self.samps if len(sm.bands)>0]
        for sm in self.samps:
            sm.wafer=self
        if position>0:
            pos=np.genfromtxt(pattern%irange[0],max_rows=2).T
            self.samps=self.transpose(pos=pos[1:])

    def corrdark(self,dark):
        for sm in self.samps:
            sm.corrdark(dark)
            
    def band_limit(self,iband=0,xmin=0,xmax=None,reset=False):
        for sm in self.samps:
            if len(sm.bands)<iband+1: continue
            bd=sm.bands[iband]
            if reset or not np.iterable(bd.sel):
                bd.sel=bd.ix>=xmin
            else:
                if xmin>0:
                    bd.sel[bd.ix<xmin]=False
            if xmax>0:
                bd.sel[bd.ix>xmax]=False
        

    def transpose(self,nbext=0,pos=[]):
        if nbext==0: nbext=len(self.samps)
        gx=[sm.bands[0].ix for sm in self.samps]
        sm=self.samps[0]
        nsamps=[]
        for j in range(len(sm.bands)-1):
            gy=[sm.bands[j].iy for sm in self.samps]
            sdata=[]
            for i in range(nbext):
                sdata.extend([gx[i],gy[i]])
            nsamps.append(Sample(None,data=sdata,maxband=nbext))
            nsamps[-1].wafer=self
            if j<len(pos):
                nsamps[-1].pos=pos[j]
        #self.make_hash()
        return nsamps
    
    def get_pos():
        return np.array([list(sm.pos) for sm in self.samps]).T
    
    def get_edge(self,pos=None,rep=1):
        if not np.iterable(pos):
            pos=self.get_pos()
        from scipy import ndimage
        xlst=np.unique(pos[0])
        ylst=np.unique(pos[1])
        out=np.array([xlst,ndimage.minimum(pos[1],pos[0],xlst)]).T
        out=np.r_[out,np.array([xlst,ndimage.maximum(pos[1],pos[0],xlst)]).T]
        out2=np.array([ndimage.minimum(pos[0],pos[1],ylst),ylst]).T
        out2=np.r_[out2,np.array([ndimage.maximum(pos[0],pos[1],ylst),ylst]).T]
        #l.scatter(*pos)
        cent=pos.mean(1)
        allang=np.r_[np.arctan2(out[:,0]-cent[0],out[:,1]-cent[1]),np.arctan2(out2[:,0]-cent[0],out2[:,1]-cent[1])]
        sorang=allang.argsort()
        dif=allang[sorang][1:]-allang[sorang][:-1]
        okang=[sorang[0]]+list(sorang[1:][dif>0])
        edgpos=np.r_[out,out2][okang]
        if rep==0:
            klist=["pos_"+str(list(p))[1:-1].replace(",","_").replace(" ","") for p in edgpos]
            return [self.samps[self.names[k]] for k in klist if k in self.names]
        return edgpos
    
    def make_hash(self,nearest=10,dshift=0.001):
        self.names={}
        allpos=[]
        for sm in self.samps:
            if len(sm.pos)<1: 
                allpos.append([0,0])
                continue
            self.names[sm.get_name()]=sm
            allpos.append(sm.pos)
        from numpy import array,arange,sqrt
        allpos=array(allpos)
        anear=[]
        for sm in self.samps:
            if len(sm.pos)<2: continue 
            dist=sqrt(((allpos-sm.pos[np.newaxis,:])**2).sum(1))
            nidx=np.argsort(dist)[1:nearest+1]
            sm.nearest.clear()
            for i in nidx:
                sm.nearest[dist[i]+dshift*i]=self.samps[i]
            anear.append(dist[nidx[0]])
        print("nearest neighbour distance %.2f +- %.2f"%(np.mean(anear),np.std(anear)))
            
    def census(self,attr,vnull=0,bmin=0,bmax=None):
        rep=[]
        for sm in self.samps:
            i=bmin
            rep.append([getattr(b,attr,vnull) for b in sm.bands[bmin:bmax]])
        return rep

    def select_fun(self,fun,iband=0,mark=None):
        rep=[]
        for sm in self.samps:
            if iband>=0:
                if fun(sm.bands[iband])==False: continue
            elif iband==-1: #any negative rejects
                for bd in sm.bands:
                    if fun(bd)==False: break
                else:
                    rep.append(sm)
                    if mark!=None:
                        sm.thick[mark]=-1
                continue
            elif iband==-2: #all negative rejects
                for bd in sm.bands:
                    if fun(bd): break
                else:
                    continue
            elif iband==-3:
                if len(sm.pos<2): continue
                if fun(sm.pos[0],sm.pos[0])==False: continue
            rep.append(sm)
            if mark!=None:
                sm.thick[mark]=-1
        return rep

    def fit(self,thguess,renow=False,bmin=0,bmax=None):
        for sm in self.samps:
            i=bmin
            for b in sm.bands[bmin:bmax]:
                if renow:
                    b.renorm(thguess)
                    b.fit(thguess,b.rat,save="b%i"%(i+1))
                else:
                    b.fit(thguess,save="b%i"%(i+1))
                i+=1
    
    def mark_empty(self,mark="proc",perc=10,min_cnt=2,iband=0,seplev=None):
        if seplev==None:
            fun2 = lambda sm:np.percentile(sm.bands[1].iy,100-perc)-np.percentile(sm.bands[1].iy,perc)
            dat2=np.array([fun2(sm) for sm in self.samps])

            prof=np.histogram(dat2,20);
            seplev=np.median(prof[1][1:][prof[0]<min_cnt])
        fun1 = lambda bd:np.percentile(bd.iy,100-perc)-np.percentile(bd.iy,perc)>seplev
        return self.select_fun(fun1,mark=mark,iband=iband)

    def dump(self,fname):
        ofile=open(fname,"w")
        for sm in self.samps:
            sm.dump_hash(ofile)
        ofile.close()

    def load(self,fname):
        import os
        if not os.path.exists(fname):
            print("file %s not found"%s)
            return
        j,k=0,0
        with open(fname) as ifile:
            for l in ifile.readlines():
                vals=l.split(':')
                if len(vals)==3:
                    j+=1
                    label=vals[0].strip()
                    if not label in self.names: continue
                    sm=self.names[label]
                    try:
                        sm.thick[vals[1].strip()]=float(vals[2])
                        k+=1
                    except:
                        print("not parsing "+l)
        return j,k
                        
    def find_problems(self,lab='first',limdif=20):
        from numpy import median
        problems=[]
        for sm in self.samps:
            if not lab in sm.thick:
                problems.append(sm)
                continue
            mvals=sm.census(lab)
            if abs(sm.thick[lab]-np.median(mvals))>limdif: 
                problems.append(sm)
        return problems

    def thick(self,patt,rep=1):
        rlist=[]
        for sm in self.samps:
            if type(patt)==list:
                for p in patt:
                    if p in sm.thick.keys():
                        rlist.append(sm.thick[p])
            else:
                for n in sm.thick.keys():
                    if n.find(patt)>=0:
                        rlist.append(sm.thick[n])
        if rep==0: return rlist
        if rep==1: return np.mean(rlist)
        return np.mean(rlist),np.std(rlist)

bbase=[]
