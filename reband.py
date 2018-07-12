import numpy as np

diel={}
minmodval=0.05
refthk=0.5#23.7 #oxide thickness on reference sample in [nm]
rescale=1
#indir="C:/Users/Admin/Documents/Lab/MOCVD/lpcvd-calib/"
indir="C:/Users/Optik/Documents/Data/Calib/"

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
        self.ix,self.iy=samp.data[loc],samp.data[loc+1]
        #self.sel=slice(10,-10)
        self.sel=self.ix>0
        self.sel[:10]=False
        self.sel[-10:]=False
        self.correct=""
        self.scale=1

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
        import profit
        epsi=[diel['ksio2'](self.ix),diel['ksi'](self.ix)]
        if not(np.iterable(dlist)): dlist=np.arange(20,300,5)
        nor=lambda d:((self.absol())[self.sel]/self.model([d],renow=True)[self.sel]).std()
        slist=[nor(d) for d in dlist]
        if rep==0: return slist
        return dlist[np.argmin(slist)]

    def renorm(self,minval=0.05,thick=None):
        if thick==None:
            thick=np.median(self.samp.thick.values())
        osel=self.sel*(self.model([thick,1,0])>minval)
        yval=self.absol()
        osel*=yval>minval
        if sum(osel)<3:
            print("too few points")
            osel=self.sel

        res=np.polyfit(self.model([thick,1,0])[osel],yval[osel],1)
        if res[0]<0:
            print("neg.scale - rejected")
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
        import profit
        if len(pars)<3:
            if renow: self.renorm(minmodval,pars[0])
            if hasattr(self,"rat"): p=[pars[0],self.rat[0],self.rat[1]]
            else: p=[pars,1,0]
        else:p=pars[:3]
        if self.samp.lay!=None and self.samp.lay.find('SiN')>=0:
            if not 'ksin' in diel:
                from scipy import interpolate as ip
                tsin=np.loadtxt(indir+"sinx_cl2.mat",unpack=True,skiprows=3)
                diel['ksin']=ip.interp1d(tsin[0],(tsin[1]+1j*tsin[2])**2)
            epsi=[diel['ksin'](self.ix),diel['ksi'](self.ix)]
        else:
            epsi=[diel['ksio2'](self.ix),diel['ksi'](self.ix)]
        return profit.plate(self.ix,epsi,[p[0]])*p[1]+p[2]

    def plot(self,amodel=False):
        from matplotlib import pyplot as pl
        pl.plot(self.ix[self.sel],self.absol()[self.sel])
        if amodel:
            thick=np.median(self.samp.thick.values())
            pl.plot(self.ix[self.sel],self.model([thick])[self.sel])

    def fit(self,inval=None,irat=[1.3,-0.15],save=None):
        from scipy import optimize as op
        resid=lambda p:sum((self.absol()[self.sel]-self.model(p)[self.sel])**2)
        zpar=op.fmin(resid,[inval,irat[0],irat[1]])
        self.qfit=resid(zpar)
        if save!=None:
            self.samp.thick[save]=zpar[0]
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

    def __init__(self,fname,laystruct="SiO2/Si",delim=None,maxband=0,data=None):
        self.bands=[]
        self.thick={}
        self.fname=fname
        if fname==None:
            assert np.iterable(data)
            self.data=data
        else:
            import os
            if not(os.path.exists(fname)):
                print("file "+fname+" not found")
                return
            self.data=np.loadtxt(self.fname,unpack=True,delimiter=delim)
            if self.data.shape[0]>self.data.shape[1]: self.data=self.data.T
        if maxband==0: maxband=len(self.data)//2
        for i in range(0,maxband*2,2):
            self.bands.append(Band(self,i))
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
                x=np.r_[0.5:4.5:0.01]#epssi[0]
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

    def plot(self,amodel=False,lims=[0,1],unit='eV'):
        for b in self.bands:
            b.plot(amodel)
        from matplotlib import pyplot as pl
        pl.ylim(*lims)
        pl.xlabel(unit)
        pl.grid()


class Wafer():

    expos=0

    def __init__(self,pattern,irange,laystruct=None,delim=None,maxband=0):
        import os
        self.samps=[Sample(pattern%i,laystruct=laystruct,delim=delim,maxband=maxband) for i in irange]
        self.samps=[sm for sm in self.samps if len(sm.bands)>0]
        for sm in self.samps:
            sm.wafer=self

    def corrdark(self,dark):
        for sm in self.samps:
            sm.corrdark(dark)

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
