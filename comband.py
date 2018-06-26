import numpy as np

def overlaps(caltab):
    cn=(np.sum([c.sum(0) for c in caltab],0)-0.9).astype(int)
    #pl.plot(enx,cn)
    istart=np.where(cn[1:]>cn[:-1])[0]
    iend=np.where(cn[1:]<cn[:-1])[0]
    return list(istart),list(iend)

def gettrans(enx,ref,dpos,smot=0.03,skiplowest=0,rep=1,scale=[],spow=1):
    '''create transformation table
    '''
    ysel=None
    if skiplowest>0:
        ysel=[ref[i]>np.percentile(ref[i],skiplowest) for i in range(len(ref))]
        ref=[ref[i][ysel[i]] for i in range(len(ref))]
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
        zub/=zn[np.newaxis,:]
        for j in range(len(caltab)):
            caltab[j][:,istart[k]:iend[k]+1]*=zub[j]
    if len(scale)==len(caltab):
        for j in range(len(caltab)):
            caltab[j]/=scale[j]
    return caltab,ysel

from scanner import spectra
enx=None#np.r_[1.2:6.4:0.01]
epssi=None#spectra.dbload("cSi_asp",enx)[1]

def sfine2(mom,leps,thk=254,emax=3.4,minor=0,dofit=4,inorm=[1,1,1],doplot=False):
    '''
    dofit==1: fit only thickness, renorm parameters fixed
    minor: minim. value of refer. spectrum to include point in chi2
    '''
    import profit
    inipars=[thk]+inorm
    enxz=enx[enx<emax]
    def zres(pars):
        model1=profit.plate(enxz,[leps,epssi[enx<emax]],pars[:1])
        dif=0
        for i in range(len(mom)):
            dsel=(mom[i]>0)*(enx<emax)
            if minor>0: dsel*=nor[i]>minor
            dif+=((mom[i][dsel]*pars[i+1]-model1[dsel[enx<emax]])**2*nor[i][dsel]).sum()*wnor[i]
        return dif
    from scipy import optimize as op
    if dofit==1: 
        fitpars=op.fmin(lambda p:zres(list(p)+inipars[1:]),inipars[:1],disp=0)
        qpars=inipars[1:]
    else: 
        fitpars=op.fmin(zres,inipars,disp=0)
        qpars=fitpars[1:]
    if doplot:
        from matplotlib import plot as pl
        model1=profit.plate(enxz,[leps,epssi[enx<emax]],fitpars[:1])
        pl.plot(enxz,model1)
        for i in range(len(mom)):
            m=mom[i]
            pl.plot(enx[(m>0)*(enx<emax)],m[(m>0)*(enx<emax)]*qpars[i],'rgm'[i])
    if dofit==1: 
        return fitpars,zres(list(fitpars)+inipars[1:])
    else: 
        return fitpars,zres(fitpars)

def comfine2(corr,irang,emax=3.4,minor=0,cord=2):
    '''dofit==1: spolecna renormalizace pro vsechny
       dofit==4: kazdy zvlast
       cord: rad korekce diel. funkce'''
    resi=[]
    enxz=enx[enx<emax]
    leps=epstio(enxz)*np.polyval(corr[:cord],enxz-enxz.mean())#(corr[0]*(enxz-enxz.mean())+corr[1])
    dofit=(len(corr)>cord) and 1 or 4
    for i in irang:
        if dofit==1: resp=cb.sfine2(momall[i],leps,results[i][0],emax=emax,inorm=list(corr[cord:]),dofit=dofit,minor=minor)
        else: resp=cb.sfine2(mom2all[i],leps,results[i][0],emax=emax,inorm=list(results[i][1:]),minor=minor)
        resi.append(resp[1])
        restemp[i][:dofit]=resp[0][:dofit]
    if loud>1:print(min(resi),max(resi))
    return sum(resi)

renor=lambda d:[(d[i]-drk[i])/nref[i] for i in range(3)]

def renorm(preps=[]):
    qcorr=np.array([ 1.40927795,  1.35854274,  1.53249378])
    if (len(preps)==0): #reload
        xpos=[np.loadtxt(indir+a)[:,0] for a in ls1 if a.find('dark')==0]
        drk=[np.loadtxt(indir+a)[:,1] for a in ls1 if a.find('dark')==0]
        ref=[np.loadtxt(indir+a)[:,1] for a in ls1 if a.find('ref')==0]
        nref=np.array(ref[:3])-np.array(drk)
        nref[nref<5]=5
        preps=[l[:8] for l in ls1 if l.find('_00')>0 and l[:2]=='ti']
        
        ratall=[renor([np.loadtxt(indir+a)[:,1] for a in ls1 if a.find(pf)==0]) for pf in preps] 
    band2,bsel2=cb.gettrans(enx,nref,xpos,smot=0.01,skiplowest=7,rep=1,scale=1/qcorr,spow=1)
    mom2all=[np.sum([rt[i][bsel2[i]].dot(band2[i])*sior for i in range(3)],axis=0) for rt in ratall]
    esel*=enx>1.35
    ok=[pl.plot(enx[esel],zoo[esel]) for zoo in mom2all[2:-1]]

    nor=[nref[i][bsel2[i]].dot(band2[i]) for i in range(len(nref))]
    wnor=[nref[i].sum()/1e4/band2[i].sum() for i in range(len(nref))]

    locfile="tio2_tl-eps.mat"
    en2,er,ei=np.genfromtxt("http://physics.muni.cz/~munz/manip/CVD/"+locfile,skip_header=3)[::-1].T
    epstio=lambda e:ip.interp1d(en2,er)(e)+1j*ip.interp1d(en2,ei)(e)

    conv1=[  2.56081863e+02,  8.21570652e-02]
    allfit=[spectra.fitting(enx[esel],zoo[esel]) for zoo in mom2all[2:-1]]
    results=np.zeros((len(momall),4))
    results[:,0]=[conv1[0]/(f[0][3]+conv1[1]) for f in allfit]
    results[:,1:]=qcorr
    restemp=np.zeros_like(results)

