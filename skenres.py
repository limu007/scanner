global slist,valrange
slist=[]
valrange=[300,550]
refsamp=None
import numpy as np
#global stacked_data
backfid='second'
prefit=False
fix_zero=True
has_bounds=True
prefun=True
import reband as rb

#preparation

def load_calib(fname='calib2.txt',iord=6,lims=[],band=0,doplot=False):
    cal1=np.loadtxt(fname)#,unpack=True)
    if len(lims)==2: sel=(cal1[2*band]>=lims[0])*(cal1[2*band]<lims[1])
    else: sel=cal1[2*band]>0
    x=cal1[2*band][sel]
    if doplot: pl.plot(x,cal1[2*band+1][sel])
    idx=np.polyfit(x,cal1[2*band+1][sel],iord)
    if doplot:
        pl.plot(x,np.polyval(idx,x))
        pl.ylim(.6,1.1)
    return idx

def recalib(waf,idx1,idx2,invert=False,ibands=[0],vlims=[0.3,1.1]):
    for j in ibands:
        x=waf.samps[0].bands[j].ix
        if np.iterable(idx1):
            corval1=np.polyval(idx1,x)
            corval1[corval1<vlims[0]]=vlims[0]
            corval1[corval1>vlims[1]]=vlims[1]
        else:
            corval1=idx1
        if invert: corval1=1/corval1
        if np.iterable(idx2):
            corval2=np.polyval(idx2,x)
            corval2[corval2<vlims[0]]=vlims[0]
            corval2[corval2>vlims[1]]=vlims[1]
        else:
            corval2=idx2
        if invert: corval2=1/corval2
        for i in range(len(waf.samps)):
            sm=waf.samps[i]
            sm.bands[j].iy/=(i*corval2+(len(waf.samps)-i)*corval1)/len(waf.samps)

#empirical fitting
def spec_fit(waf,prange=[.2,2],lbin=8,loud=0,istart=50,uband=0,esel=[],mode='inve'):
    bd=waf.samps[istart].bands[uband]
    zsel=bd.sel.copy()
    if len(esel)>1:
        zsel[bd.ix<esel[0]]=False
        zsel[bd.ix>esel[1]]=False
    x=bd.ix[zsel].copy()
    x-=x.mean()
    import spectra
    of,chi=spectra.fitting(x,bd.absol()[zsel],lbin=lbin,prange=prange,loud=abs(loud),fit_mode=0,refr_mode=mode)#"simpl")
    apars,achi=[],[]
    if loud<0: return of,chi
    for s in waf.samps:
        cof,chi=spectra.fitting(x,s.bands[uband].absol()[zsel],lbin=lbin,p0=of,prange=prange,loud=0,fit_mode=0,refr_mode=mode)
        apars.append(cof)
        achi.append(chi)
    return np.array(apars).T,np.array(achi)

def fitme(sm,lab="first",inlab=None,loud=0):
    if inlab==None: inlab=lab
    cvals=sm.census(inlab,valid=valrange)
    if len(cvals)>0:
        thk=np.median(cvals)
    else:
        thk=(valrange[0]+valrange[1])/2.#waf3.thick(inlab)
    if thk<valrange[0]: thk=valrange[0]
    if thk>valrange[1]: thk=valrange[1]
    if loud>0:
        print("thick %.1f nm"%thk)
    if has_bounds:
        bounds=[valrange]
        for i in range(len(sm.bands)):
            if fix_zero: bounds+=[[rb.min_norm,rb.max_norm]]
            else: bounds+=[[rb.min_norm,rb.max_norm],[-rb.max_bias,rb.max_bias]]
    else:
        bounds=None
    res=sm.fit(thk,save=lab,refer=refsamp,prefit=prefit,prefun=prefun,fix_zero=fix_zero,constr=bounds)
    #if fix_zero: sm.renorm()
    return res

def fitlist(lst):
    return [fitme(slist[l]) for l in lst]

def save_id(s,vals,fid):
    try:
        sm=slist[s]
    except:
        sel=[sm for sm in slist if s==sm.get_name()]
        if len(sel)==0: 
            print(str(s)+" not found")
            return
        sm=sel[0]
    sm.thick[fid]=vals[0]
    for k in range(len(sm.bands)):

        if len(vals)>2*k+2:
            sm.bands[k].rat=vals[2*k+1:2*k+3]
        rat=sm.bands[k].rat
        if refsamp!=None and len(refsamp.bands)>k:
            sm.bands[k].qfit=slist[s].bands[k].fit(None,refer=refsamp.bands[k])([vals[0],rat[0],rat[1]]) #quality of fit
            sm.bands[k].qfit=slist[s].bands[k].fit(None)([vals[0],rat[0],rat[1]])
        #slist[s].bands[k].qfit=slist[s].bands[k].fit(None)(val) 
    return sm

def inpool(cores_used=5,fid='first',offs=0,missing=False,gsiz=10,rewrite=False):

    if missing: olist=[i for i in range(len(slist)) if not fid in slist[i].thick]
    else: olist=np.r_[offs:len(slist)]

    maxi=(len(olist)-offs)//gsiz-1
    wlist=[olist[i*gsiz:i*gsiz+gsiz] for i in range(maxi)]
    if (gsiz*maxi)<len(olist):
        wlist.append(olist[gsiz*maxi:])
    from multiprocessing import Pool
    p=Pool(processes = cores_used)
    stacked_data = p.map(fitlist,wlist)

    for j in range(len(wlist)):
        for i in range(len(wlist[j])):
            save_id(wlist[j][i],stacked_data[j][i],backfid if (not rewrite and fid in slist[s].thick) else fid)
    return stacked_data

def clean_out(waf,th2,sd2,nsig=2,lab='both'):
    k=0
    for i in range(len(waf.samps)):
        s=waf.samps[i]
        if lab in s.thick:
            if s.thick[lab]<th2-nsig*sd2 or s.thick['both']>th2+nsig*sd2:
                del s.thick[lab]
                k+=1
    return k

def iterate(glist,lab="both",doplot=True,margin=50,gsiz=8,maxiz=0,thread=True):
    #th2=waf3.thick(lab)
    #global valrange
    #valrange=[th2-margin,th2+margin]
    '''selects samples with some results in neighborhood
    '''
    olist=[s for s in glist if len(s.census(lab,valid=valrange))>0 and not lab in s.thick]
    if maxiz>0: olist=olist[:maxiz]
    #if len(slist)==0:
    #    return th2,0
    if doplot:
        from matplotlib import pyplot as pl
        rep2y=np.array([list(sm.pos) for sm in olist])
        #pl.scatter(rep2y[:,0],rep2y[:,1],marker="x")
    if gsiz<0: return olist
    if thread:
        global slist
        slist=olist
        sethk=inpool(fid=lab,missing=False,rewrite=True,gsiz=gsiz)[0]
    else:
        sethk=[fitme(s)[0] for s in olist]
        if len(sethk)>0:
            for i in range(len(olist)):
                olist[i].thick[lab]=sethk[i]
            return len(olist),min(sethk),max(sethk)
        else:
            return len(olist),len(sethk)
    return len(olist)#,min(set3),max(set3)
#stacked_data=inpool(slist,fid="both",missing=True,rewrite=True)

def process_rest(glist,gsiz=5,lab="both",cores_used=6,maxiter=15,minsamps=5,doplot=False,thread=True):
    '''uses ssamps
    '''
    global slist
    for k in range(maxiter):
        slist=iterate(glist,lab=lab,gsiz=-5,thread=False,doplot=doplot) #ziskani seznamu
        if len(slist)<=minsamps or not thread:
            wlist=[fitme(s,lab=lab)[0] for s in slist]
            if len(slist)<=minsamps: break
            continue
        olist=np.r_[:len(slist)]
        maxi=len(olist)//gsiz-1
        wlist=[olist[i*gsiz:i*gsiz+gsiz] for i in range(maxi)]
        if (gsiz*maxi)<len(olist):
            wlist.append(olist[gsiz*maxi:])
        from multiprocessing import Pool
        print(len(slist))
        p=Pool(processes = cores_used)
        stacked_data = p.map(fitlist,wlist)
        for j in range(0,len(wlist)):
        #set3=Out[45]
            for i in range(len(wlist[j])):
                s=wlist[j][i]
                save_id(s,stacked_data[j][i],lab)
    #iterate(gsiz=5,thread=False)
    
def wafplot(waf,col=2,lims=[450,485],lab='both',fname=None):
    '''
    columns:
    2: thickness
    3: quality band 0
    4: quality band 1 # if more than 1 band
    5,6: linear correction band 0
    7,8: linear correction band 1 # if more than 1 band
    '''
    from matplotlib import pyplot as pl
    if len(waf.samps[0].bands)>1:
        rep2=np.array([list(sm.pos)+[sm.thick[lab],sm.bands[0].qfit,sm.bands[1].qfit]+list(sm.bands[0].rat)+list(sm.bands[1].rat) for sm in waf.samps if lab in sm.thick])
    else:
        rep2=np.array([list(sm.pos)+[sm.thick[lab],sm.bands[0].qfit]+list(sm.bands[0].rat) for sm in waf.samps if lab in sm.thick])
    rep2b=np.array([list(sm.pos)+[0] for sm  in waf.samps if not lab in sm.thick])
    pl.scatter(rep2b[:,0],rep2b[:,1],marker="x")
    if len(lims)>0:
        pl.scatter(rep2[:,0],rep2[:,1],c=rep2[:,col],vmin=lims[0],vmax=lims[1])
    else:
        pl.scatter(rep2[:,0],rep2[:,1],c=rep2[:,col])
    pl.colorbar()
    if fname!=None:pl.savefig(fname)
