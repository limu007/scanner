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

def fitme(sm,lab="first",inlab=None):
    if inlab==None: inlab=lab
    cvals=sm.census(inlab,valid=valrange)
    if len(cvals)>0:
        thk=np.median(cvals)
    else:
        thk=(valrange[0]+valrange[1])/2.#waf3.thick(inlab)
    if thk<valrange[0]: thk=valrange[0]
    if thk>valrange[1]: thk=valrange[1]
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
    4: quality band 1
    5,6: linear correction band 0
    7,8: linear correction band 1
    '''
    from matplotlib import pyplot as pl
    rep2=np.array([list(sm.pos)+[sm.thick[lab],sm.bands[0].qfit,sm.bands[1].qfit]+list(sm.bands[0].rat)+list(sm.bands[1].rat) for sm in waf.samps if lab in sm.thick])
    rep2b=np.array([list(sm.pos)+[0] for sm  in waf.samps if not lab in sm.thick])
    pl.scatter(rep2b[:,0],rep2b[:,1],marker="x")
    if len(lims)>0:
        pl.scatter(rep2[:,0],rep2[:,1],c=rep2[:,col],vmin=lims[0],vmax=lims[1])
    else:
        pl.scatter(rep2[:,0],rep2[:,1],c=rep2[:,col])
    pl.colorbar()
    if fname!=None:pl.savefig(fname)
