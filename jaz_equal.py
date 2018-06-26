renor=lambda d:[(d[i]-drk[i])/nref[i] for i in range(3)]

preps=[]
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


tiocorr=np.array([-0.06273044,  1.09829889])
conv1=[  2.56081863e+02,  8.21570652e-02]
allfit=[spectra.fitting(enx[esel],zoo[esel]) for zoo in mom2all[2:-1]]
results=np.zeros((len(momall),4))
results[:,0]=[conv1[0]/(f[0][3]+conv1[1]) for f in allfit]
results[:,1:]=qcorr
restemp=np.zeros_like(results)
loud=1
emax=3.1
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
results=restemp.copy()
corr33=op.fmin(comfine,list(tiocorr)+list(qcorr),(range(len(mom2all)),emax,0.1,3))
corr33
