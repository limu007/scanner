import numpy as np

refer=None
diel={}
minmodval=0.05
min_normstep,min_method=0.001,'L-BFGS-B'#'SLSQP'
min_norm,max_norm=0.5,1.7 #norm correction
max_bias=0.45
refthk=0.5#23.7 #oxide thickness on reference sample in [nm]
rescale=1
#indir="C:/Users/Admin/Documents/Lab/MOCVD/lpcvd-calib/"
indir="C:/Users/Optik/Documents/Data/Calib/"
sinname="sin_nk_test.mat"#"sinx_cl2.mat"

def transpose(self,nbext=0,pos=[]):
    '''switching columns/files in data: channels vs. measurement points (Samples)'''
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

def testgap(cens,limfrac=0.1,mindif=10):
    if len(cens)<2: return
    cens=np.array(cens)
    mval,sval=np.mean(cens),np.std(cens)
    if sum(abs(cens-mval)<0.5*sval)/len(cens)<limfrac: #has gap
        co1=cens<mval-.5*sval
        mval1,sval1=np.median(cens[co1]),np.std(cens[co1])
        co2=cens>mval+.5*sval
        mval2,sval2=np.median(cens[co2]),np.std(cens[co2])
        if (mval2-mval1)<mindif: return mval,sval
        return mval1,sval1/sval,mval2,sval2/sval
    return mval,sval

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
        #from scanner import profit
        #epsi=[diel['ksio2'](self.ix),diel['ksi'](self.ix)]
        if not(np.iterable(dlist)): dlist=np.arange(20,300,5)
        nor=lambda d:((self.absol())[self.sel]/self.model([d],renow=True)[self.sel]).std()
        slist=[nor(d) for d in dlist]
        if rep==0: return slist
        return dlist[np.argmin(slist)]

    def renorm(self,minval=0.05,thick=None,ngrp=15.,minstep=0,loud=0):
        if minstep==0:
            minstep=min_normstep
        if thick==None:
            if len(self.samp.thick)>0:
                thick=self.samp.get_thick(unc=False)#np.median(list(self.samp.thick.values()))
            else: return
        modval=self.model([thick,1,0])
        modval[np.isnan(modval)]=0 # proc tu jsou vadne hodnoty?
        osel=self.sel*(modval>minval)
        yval=self.absol()
        osel*=yval>minval
        if sum(osel)<3:
            print("too few points")
            osel=self.sel.copy()

        #robust init
        from scipy import ndimage as nd
        mmin=modval[osel].min()
        mstep=(modval[osel].max()-mmin)/ngrp
        if mstep<minstep/3.:
            ngrp=int((modval[osel].max()-mmin)/minstep)
            if ngrp>2:
                mstep=(modval[osel].max()-mmin)/ngrp
        if mstep<minstep/3.:
            if loud>0:
                print("not variable enough [%i / %.2g]"%(ngrp,mstep))
            self.rat=[1,0.01]
            return
        modlab=((modval[osel]-mmin)/mstep).astype(int)
        #meds=nd.median(yval[osel],modlab,range(max(modlab)-1))
        #meds=nd.mean(yval[osel],modlab,range(max(modlab)-1))
        meds=np.array([np.mean(yval[osel][modlab==i]) for i in range(max(modlab)-1)])
        if loud>1:
            print(meds)
            cnt=np.array([sum(modlab==i) for i in range(max(modlab)-1)])
            print(cnt)
            print(modlab.min(),modlab.max(),sum(modlab==modlab.max()))
            print(nd.standard_deviation(yval[osel],modlab,range(max(modlab)-1)))
        #first a robust fit
        vpos=mmin+mstep*np.r_[0.5:ngrp]
        vpos=vpos[:len(meds)]
        res0=np.polyfit(vpos,meds,1)
        if loud>0:
            print(res0)
        dif=abs(modval[osel]*res0[0]+res0[1]-yval[osel])
        lim=np.percentile(dif,90)
        dif[np.isnan(dif)]=2*lim 
        isel=np.r_[:len(osel)][osel][dif>lim] #indices where residual is too big
        osel[isel]=False
        res=np.polyfit(modval[osel],yval[osel],1)
        if loud>0:
            chi2=np.sum((np.polyval(res,modval[osel])-yval[osel])**2)
            chi2out=np.sum((np.polyval(res,modval[isel])-yval[isel])**2)
            print(res,chi2,chi2out)
        if (res[0]<min_norm):
            res[0]=min_norm
        if (res[0]>max_norm):
            #print("neg.scale - rejected")
            res[0]=max_norm
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

    def model(self,pars,renow=False,prefun=False):
        '''currently supports single layer model
           dielectric as SiN, SiO2 and cauchy models
           substrate is Si
           TODO: more variability!

           prefun: returns reflectivity with precalculated dielectric functions
        '''
        from . import profit
        if len(pars)<3:
            if renow: self.renorm(minmodval,thick=pars[0])
            if hasattr(self,"rat"): p=[pars[0],self.rat[0],self.rat[1]]
            else: p=[pars,1,0]
        else:p=pars[:3]
        substr=diel['ksi'](self.ix)
        if self.samp.lay!=None:
            tlay=self.samp.lay.split('/')[0]
            if tlay=='SiN':#self.samp.lay.find('SiN')>=0:
                if not 'ksin' in diel:
                    from scipy import interpolate as ip
                    import os
                    if os.path.exists(indir+sinname):
                        tsin=np.loadtxt(indir+sinname,unpack=True,skiprows=3)
                    else:
                        print("cannot find SiN data: "+indir+sinname)
                    diel['ksin']=ip.interp1d(tsin[0],(tsin[1]+1j*tsin[2])**2)
                epsi=[diel['ksin'](self.ix),substr]
            elif tlay=='cau': #Caucjy profile
                if len(pars)>3:
                    if not 'cau' in diel: diel['cau']=lambda x:np.polyval(pars[3:][::-1],x**2)**2
                    epsi=[np.polyval(pars[3:][::-1],self.ix**2)**2,diel['ksi'](self.ix)]
                elif 'cau' in diel:
                    epsi=[diel['cau'](self.ix),substr]
            elif tlay in diel: #user defined dielectric function
                epsi=[diel[tlay](self.ix),substr]
            else:#expecting SiO2
                epsi=[diel['ksio2'](self.ix),substr]
        else: #by default
            epsi=[diel['ksio2'](self.ix),substr]
        if prefun:
            if hasattr(self,"rat") and len(pars)<2: #fixed renormalization
                return lambda p:profit.plate(self.ix,epsi,[p[0]])*self.rat[0]+self.rat[1]
            return lambda p:profit.plate(self.ix,epsi,[p[0]])*p[1]+p[2]
        return profit.plate(self.ix,epsi,[p[0]])*p[1]+p[2]

    def plot(self,amodel=False,match=True,ax=None):
        if ax==None:
            from matplotlib import pyplot as pl
        else:
            pl=ax
        if match:
            dat=pl.plot(self.ix[self.sel],(self.absol()[self.sel]-self.rat[1])/self.rat[0])[0]
        else:
            dat=pl.plot(self.ix[self.sel],self.absol()[self.sel])[0]
        if amodel and len(self.samp.thick)>0:
            thick=self.samp.get_thick(unc=False)
            if match:
                ## proc, kdyz to umi spocitat model uvnitr? 
                ## chceme model oprosteny od renormalizace..aby navazovaly sousedni intervaly(?)
                modl=pl.plot(self.ix[self.sel],(self.model([thick])[self.sel]-self.rat[1])/self.rat[0])[0]
                #pl.plot(self.ix[self.sel],self.model([thick,1,0])[self.sel])
            else:
                modl=pl.plot(self.ix[self.sel],self.model([thick])[self.sel])[0]
            return [dat,modl]
        return [dat]

    def fit(self,inval=None,irat=[1.3,-0.15],save=None,refer=None,prefun=False,constr=None):
        '''
        irat: expected values for renormalization
        save: label for storing fit results (in self.sample)
        '''
        from scipy import optimize as op
        if prefun: # save
            iabs=self.absol()[self.sel]
            imod=self.model([0]+irat,prefun=True)
            if refer!=None:
                resid=lambda p:sum((iabs-imod(p)[self.sel])**2/refer.iy[self.sel]**2)
            else:
                resid=lambda p:sum((iabs-imod(p)[self.sel])**2)
            return resid
        if refer!=None:
            resid=lambda p:sum((self.absol()[self.sel]-self.model(p)[self.sel])**2/refer.iy[self.sel]**2)
        else:
            resid=lambda p:sum((self.absol()[self.sel]-self.model(p)[self.sel])**2)
        if inval==None: return resid

        from numpy import iterable
        if iterable(inval):
            initv=[inval[0],irat[0],irat[1]]+list(inval[1:])
        else:
            initv=[inval,irat[0],irat[1]]
        if iterable(constr):
            res=op.minimize(resid, initv, method=min_method, bounds=constr)
            zpar=res.x
            if hasattr(res,'hess_inv'):
                if hasattr(res.hess_inv,'todense'):
                    cov=res.hess_inv.todense()
                    self.err=np.sqrt(cov.diagonal())
                    self.cor=cov/self.err[:,np.newaxis]/self.err[np.newaxis,:]
        else:
            zpar=op.fmin(resid,initv,full_output=False,disp=False)
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
    data=None
    fname=""
    lay=None
    wafer=None
    pos=[]
    inval_thick=set()#[]
    normtab={}

    def __init__(self,fname,laystruct="SiO2/Si",delim=None,maxband=0,data=None,headerow=0,bord=10):
        self.bands=[]
        self.thick={}
        self.thickerr={}
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
            if maxband==-1:
                maxband=-len(self.data)+1 # columns = [energy val_pos_1 val_pos_2]
            for i in range(-1,maxband,-1):
                self.bands.append(Band(self,i,bord=bord))
        else:
            for i in range(0,maxband*2,2): #columns = [energy value energy_2 value_2 ...]
                self.bands.append(Band(self,i,bord=bord))
        #self.norm=lambda x:1
        if laystruct!=None: self.lay=laystruct

    def calib(self,rc=None):
        from scanner import spectra
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
        if not 'ksin' in diel:
            if hasattr(rc,'ref_epsweb'):
                try:
                    tsin=np.loadtxt(rc.ref_epsweb+rc.ref_epsfile['SiN'],skiprows=3,unpack=True)
                except:
                    print('cannot fetch from ',rc.ref_epsweb)
                else:
                    diel['ksin']=ip.interp1d(tsin[0],(tsin[1]+1j*tsin[2])**2)
        if not 'ksio2' in diel:
            if not os.path.exists(indir+"sio2_palik_g.mat"):
                x=np.r_[0.5:6.5:0.01]#epssi[0]
                tsio2=[x,np.polyval(spectra.cau_sio2,x)]
            else:
                tsio2=np.loadtxt(indir+"sio2_palik_g.mat",unpack=True,skiprows=3)
            diel['ksio2']=ip.interp1d(tsio2[0],tsio2[1]**2)
        if self.wafer!=None: self.wafer.status['calibrated']=True

    def norm(self,px):
        from scanner import profit
        bid='b%.1f-%.1f'%(px[0],px[-1])
        if not bid in self.normtab:
            if refthk<1.: self.normtab[bid] = profit.reflect(diel['ksi'](px))
            else: self.normtab[bid] = profit.plate(px,[diel['ksio2'](px),diel['ksi'](px)],[refthk])
        return self.normtab[bid]
        
    def corrdark(self,dark):
        for i in range(len(self.bands)):
            if i>=len(dark.bands): break
            self.bands[i].corrdark(dark.bands[i].iy)
    
    def update_nearest(self,maxnbrh=10,maxdist=25):
        import numpy as np
        cpos=np.array(self.pos)
        smsel=[s for s in self.wafer.samps if max(abs(cpos-s.pos))<maxdist]
        dist=[np.sqrt(sum((cpos-s.pos)**2)) for s in smsel]
        zlist=np.argsort(dist)
        if maxnbrh>=len(smsel): maxnbrh=len(smsel)-1
        ddist=0.001
        self.nearest=nnear=dict([(dist[j]+ddist*j,smsel[j]) for j in zlist[1:maxnbrh+1]])
        return np.mean(dist[1:maxnbrh+1])

    def get_predict(self,pord=2,trange=[]):
        '''interpolate already calculated valued to get local prediction
        '''
        pos=np.array(self.pos)
        tab=[]
        for s in self.nearest.values():
            thk=s.get_thick()
            if thk==None: continue
            tab.append(np.r_[pos-s.pos,thk.n])
        predtab=np.array(tab).T
        if len(trange)>1:
            sel=(predtab[2]>=trange[0])*(predtab[2]<=trange[1])
            predtab=predtab[:,sel]
        px,py,predval=predtab
        if len(px)<4: return #cannot fit anythin`
        if (pord==2) and len(px)>6: modA=np.array([np.ones_like(px),px,py,px**2,py**2,px*py]) #quadratic
        else:modA=np.array([np.ones_like(px),px,py])
        matH=modA@modA.T
        predic=np.linalg.inv(matH)@modA@predval
        chi2=((predval-modA.T@predic)**2).sum()
        return predic[0],chi2

    def get_thick(self,select=None,unc=True):
        if len(self.thick)==0: return None
        vlist=[]
        if select!=None:
            klist=[k for k,v in self.thick.items() if k.find(select)>=0]
        elif len(self.inval_thick)>0:
            klist=[k for k,v in self.thick.items() if k not in self.inval_thick]
        if len(vlist)==0:
            klist=list(self.thick.keys())
        vlist=np.array([self.thick[k] for k in klist if self.thick[k]>0])
        if len(vlist)==0: return None
        elist=np.array([self.thickerr.get(k,0) for k in klist if self.thick[k]>0])
        if np.all(elist==0): elist[:]=1
        else: elist[elist==0]=elist.max()
        if len(vlist)==1:
            vcom,ecom=vlist[0],elist[0]
        else:
            ecom=1/np.sqrt((1/elist**2).sum())
            vcom=((vlist/elist**2).sum()*ecom**2)
        if not unc: return vcom
        import uncertainties as uc
        return uc.ufloat(vcom,ecom)
        #return np.median(vlist)

    def get_qfit(self,thk=None,save=True):
        res=[]
        if thk==None:
            thk=self.get_thick(unc=False)
            if thk==None: return
        for bd in self.bands:
            if not hasattr(bd,'rat'):
                bd.renorm(thick=thk)
            bresid=bd.fit(None)
            res.append(bresid([thk]+list(bd.rat)))
            if save: bd.qfit=res[-1]
        return res

    def get_der2(self,dis=0.01,resf=None,arrf=None,rep=0,lab=None):
        if resf==None: resf=self.fit(prefun=True)
        if not np.iterable(arrf):
            arrf=self.fit(prefun=False)
            if lab in self.thick: arrf[0]=self.thick[lab]
            else: arrf[0]=self.get_thick(None,False)
            arrf=arrf.astype(float)
        avec=np.eye(len(arrf))
        
        chimin=resf(arrf)
        upvec=[resf(arrf*(1+dis*av)) for av in avec]
        dnvec=[resf(arrf*(1-dis*av)) for av in avec]
        dder=(np.array(upvec)+dnvec-2*chimin)/dis**2/arrf.astype(float)**2
        if rep==2: return np.sqrt(1/dder),resf,arrf
        if rep==1: return np.sqrt(1/dder),arrf
        return np.sqrt(1/dder) #sigma v rezu
    

    def model_cov(self,msiz=40,dis=0.01,minder=0,lab=None):
        qder,resf,arrf=get_der2(self,rep=2,lab=lab)
        #qder=(np.array(upvec)+dnvec-2*chim)/dis**2
        sel=qder>minder
        amat=np.array([sel*np.random.normal(size=len(sel)) for i in range(msiz)])
        #mozna upravit velikost kroku v tom kterem smeru
        #Hmat=amat.T@amat
        #covX=np.linalg.inv(Hmat[sel][:,sel])
        arrf=arrf.astype(float)
        amat=amat*qder*dis
        #print(abs(amat).max(0))
        zlst=np.r_[:len(sel)][sel]
        out=[]
        for j in zlst:
            out.extend([(i,j) for i in zlst[zlst>j]])
        extmat=[amat[:,ia]*amat[:,ib] for ia,ib in out]
        fullmat=np.r_[[np.ones(len(amat))],amat[:,sel].T,amat[:,sel].T**2,extmat] 

        achis=[resf(arrf+amat[i]) for i in range(len(amat))]
        #pl.plot(np.sqrt((amat**2).sum(1)),achis-chim,'d')
        fhesX=fullmat@fullmat.T
        fcovX=np.linalg.inv(fhesX)
        pars=fcovX@fullmat@achis
        sig=np.sqrt(pars[0]/(sum(sum([b.sel for b in self.bands]))-len(sel)))
        
        hesF=np.eye(len(sel))*1e-4
        for i in zlst:
            hesF[i,i]=pars[1+sum(sel)+i]
        for k in range(len(out)):
            i,j=out[k]
            hesF[i,j]=pars[1+sum(sel)*2+k]/2.
            hesF[j,i]=pars[1+sum(sel)*2+k]/2.
        return np.linalg.inv(hesF)*sig**2
    
    def renorm(self,thick=None):
        for bd in self.bands:
            bd.renorm(thick=thick)
            
    def get_cov(self,pars,refer=None,fix_zero=False):
        from scipy import optimize as op
        if refer!=None:
            sigma=np.concatenate([refer.bands[i].iy[self.bands[i].sel] for i in range(len(self.bands))])
        else:
            sigma=None
        if fix_zero:
            resid=lambda p:np.concatenate([self.absol()[bd.sel]-self.model([p[0],p[i+1],0])[bd.sel] for i,bd in enumerate(self.bands)])
        else:
            resid=lambda p:np.concatenate([self.absol()[bd.sel]-self.model([p[0],p[2*i+1],p[2*i+2]])[bd.sel] for i,bd in enumerate(self.bands)])
        return resid,sigma
    

    def fit(self,inval=None,save=None,refer=None,prefit=False,prefun=False,fix_zero=False,constr=None):
        '''
        prefun: returns model function, not final optimization result
        fix_zero: only slope is fitted
        '''
        from scipy import optimize as op
        from numpy import iterable
        if refer!=None:
            bresid=[self.bands[i].fit(None,refer=refer.bands[i],prefun=prefun) for i in range(len(self.bands))]
        else:
            bresid=[self.bands[i].fit(None,prefun=prefun) for i in range(len(self.bands))]
        if fix_zero:
            resid=lambda p:sum([bresid[i]([p[0],p[i+1],0]) for i in range(len(self.bands))])
        else:
            resid=lambda p:sum([bresid[i]([p[0],p[2*i+1],p[2*i+2]]) for i in range(len(self.bands))]) #sum of chi2 in all channels
        if prefun and inval==None:
            return resid
        if prefit: # fit indiviudal channels before combining all
            if refer!=None:
                bresult=[b.fit(inval,refer=refer.bands[i]) for i,b in enumerate(self.bands)]
            else:
                bresult=[b.fit(inval) for b in self.bands]
            inthk=np.median([r[0] for r in bresult]) #thickness
            inarr=np.concatenate([[inthk]]+[r[1:3] for r in bresult])
        else:
            inarr=np.concatenate([[inval]]+[b.rat for b in self.bands])
        if fix_zero:
            inarr=np.concatenate([inarr[:1],inarr[1::2]])

        if inval==None: return inarr
        # we need parameter constraints
        if iterable(constr):
            res=op.minimize(resid, inarr, method=min_method, bounds=constr)
            zpar=res.x
            #if save: return res
            if hasattr(res,'hess_inv'):
                if hasattr(res.hess_inv,'todense'):
                    cov=res.hess_inv.todense()
                    self.err=np.sqrt(cov.diagonal())
                    self.cor=cov/self.err[:,np.newaxis]/self.err[np.newaxis,:]

        else:
            zpar=op.fmin(resid,inarr,full_output=False,disp=False)
        if np.isnan(zpar[0]):
            return zpar
        if save!=None:
                self.thick[save]=zpar[0]
                self.chi2[save]=resid(zpar)
        for i,b in enumerate(self.bands):
            p=zpar
            if fix_zero:
                b.qfit=bresid[i]([p[0],p[i+1],0])
            else:
                b.qfit=bresid[i]([p[0],p[2*i+1],p[2*i+2]])
            if save!=None:
                if fix_zero:
                    b.rat[0]=p[i+1]
                else:
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

    def todo(self,lab='both'):
        return [s for s in self.nearest.values() if lab not in s.thick]

    def census(self,lab='first',outliers=20,valid=[]):
        '''colects data from (already analyzed) nearby samples
        '''
        from numpy import percentile
        rvals=[sq.thick[lab] for sq in self.nearest.values() if hasattr(sq,"thick") and lab in sq.thick and sq.thick[lab]>0]
        if len(valid)>1:
            rvals=[t for t in rvals if t>=valid[0] and t<=valid[1]]
        if len(rvals)>2 and outliers>0:
            zmin,zmax=percentile(rvals,[outliers,100-outliers])
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

    def plot(self,amodel=False,lims=[0,1],unit='eV',match=True,ax=None):
        plst=[]
        for b in self.bands:
            plst+=b.plot(amodel,match=match,ax=ax)
        from matplotlib import pyplot as pl
        pl.ylim(*lims)
        pl.xlabel(unit)
        pl.grid()
        return plst

#posfun=lambda ix,iy:np.any([(edgepos==[ix,iy]).sum(1)==2])

class Wafer():

    expos=0
    status={}

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

    def get_pos(self):
        return np.array([list(sm.pos) for sm in self.samps]).T

    def get_nearest(self,x,y,cnt=1,rep=1):
        cent=np.array([x,y]).reshape(2,1)
        dist2=((self.get_pos()-cent)**2).sum(0)
        if len(dist2)==0: return []
        smbest=self.samps[np.argmin(dist2)]
        if cnt==1: 
            if rep==2: return [[smbest,np.sqrt(dist2.min())]]
            return [smbest]
        smsom=[smbest]+[s for s in smbest.nearest.values()]
        smset=set(smsom)
        for s in smsom[1:]:
            for s2 in s.nearest.values():
                smset.add(s2)
        smsom=list(smset)
        dist2=((np.array([list(sm.pos) for sm in smsom]).T-cent)**2).sum(0)
        print(len(dist2),np.sqrt(max(dist2)))
        smids=np.argsort(dist2)[:cnt]
        if rep==2: return [(smsom[i],np.sqrt(dist2[i])) for i in smids]
        return [smsom[i] for i in smids]

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

    def mark_empty(self,mark="proc",perc=10,min_cnt=2,iband=0,rep=0,ndiv=20):
        fun2 = lambda sm:np.percentile(sm.bands[iband].iy[sm.bands[iband].sel],100-perc)-np.percentile(sm.bands[iband].iy[sm.bands[iband].sel],perc)
        dat2=np.array([fun2(sm) for sm in self.samps])
        drange=[np.percentile(dat2,10),np.percentile(dat2,90)]
        ddist=drange[1]-drange[0]
        drange[0]-=ddist
        drange[1]+=ddist
        prof=np.histogram(dat2,np.r_[drange[0]:drange[1]:1j*ndiv])
        if rep==1: return prof
        qsel=prof[0]>min_cnt
        zmin,zmax=np.r_[:len(qsel)][qsel][[0,-1]]
        seplev=np.median(prof[1][zmin:zmax+1][prof[0][zmin:zmax+1]<min_cnt])
        if rep>1: print(seplev)
        fun1 = lambda bd:np.percentile(bd.iy[bd.sel],100-perc)-np.percentile(bd.iy[bd.sel],perc)>seplev
        return self.select_fun(fun1,mark=mark,iband=iband)

    def dump(self,fname):
        ofile=open(fname,"w")
        for sm in self.samps:
            sm.dump_hash(ofile)
        ofile.close()

    def load(self,fname):
        import os
        if not os.path.exists(fname):
            print("file %s not found"%fname)
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
                        sm.thick[vals[1 ].strip()]=float(vals[2])
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
