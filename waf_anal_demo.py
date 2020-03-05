from matplotlib import pyplot as pl
import sys; sys.path.append("C:/Users/Admin/Code/")
import spectra,profit
import numpy as np
indir="C:/Users/Admin/Cloud/Lab/Spectra/Elias/preetch/"
indir="/home/limu/Public/Cloud/Lab/Spectra/Elias/preetch/"
#ifile=indir+"sio2_400nm_wf3_map_3.txt.gz"
from scanner import reband
waf=reband.Wafer(indir+"sio2_400nm_wf4_map_%i.txt",[1,2],laystruct="SiO2/Si",maxband=-374,headerow=2,position=2)
#from scanner.reband 
#from reband import transpose
#waf3.samps=transpose(waf3)
def init_lims(waf):
    sm=waf.samps[0]
    sm.calib()
    waf.make_hash()
    waf.band_limit(0,1.65,2.55)
    waf.band_limit(1,1.33,1.68)
    sm.inval_thick.add("first")#,"second","proc"]
init_lims(waf)

ssamps=waf.mark_empty(mark="first",min_cnt=2,iband=0,seplev=.5)
#ssamps=waf3.mark_empty(mark="first")
refsamp=reband.Sample(None,data=np.loadtxt("refer_waf2.dat"))
sm=ssamps[20]
#sm.calib()
th=sm.bands[0].guess(np.r_[300:500:10])
sm.fit(th,save="both",prefit=True)#,refer=refsamp)
#sm.plot(amodel=True)

import skenres
skenres.valrange=[300,550]
skenres.refsamp=refsamp

slist0=sm.nearest.values()
slist=[s for s in slist0 if len(s.census("both"))>0]
skenres.slist=slist
set1=skenres.inpool(fid="both",missing=True,rewrite=True)[0]

skenres.clean_out(waf,(skenres.valrange[0]+skenres.valrange[1])/2,50,2)

skenres.process_rest(ssamps,maxiter=5,doplot=False,minsamps=2,lab='both')
skenres.wafplot(waf)
