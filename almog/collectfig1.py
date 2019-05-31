import matplotlib
matplotlib.use('Agg')
import numpy
from pylab import *
import mytools
import pickle
import time
import sys
import random
import scipy.io

import mutation_stuff
MT = mutation_stuff.getMT()
geneNames = mutation_stuff.getgenenames()
defVals = mutation_stuff.getdefvals()

unpicklefile = open('scalings.sav', 'r')
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()
theseCoeffsAll = unpickledlist[0]
theseMutValsAll = unpickledlist[2]

styles = ['b-','b-','b-','b-','b-','b-','b-','b-','b-','b-']
#cols = ['#666666','#012345','#aa00aa','#bbaa00','#ee6600','#ff0000', '#00aaaa','#772277','#00cc00']
cols = ['#444444','#012345','#aa00aa','#bbaa00','#ee6600','#ff0000', '#009999','#772277','#00cc00']
col_control = '#2222ff'
coeffCoeffs = [[0.25,0],[0.125,0],[0.5,0],[0.5,1.0/3],[0.5,2.0/3],[0.5,1.0],[-0.25,0],[-0.125,0],[-0.5,0]]
lw = 1.5
fs = 10
ispDef = 0 # Consider a local maximum above -50mV a spike only if after last spike the membrane potential came below -50mV
#ispDef = 1 # Consider every local maximum above -50mV a spike
##Is = [0.65+0.025*x for x in range(0,11)]
#Is = [0.7+0.0125*x for x in range(0,17)]
Is = [0.7+0.005*x for x in range(0,41)]




unpicklefile = open('ifcurves_0_0_0.sav', 'r')
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()
spTimesThisMutVal = unpickledlist[1+ispDef]
nSpikes_control = [sum([1 for x in spTimesThisMutVal[5][j] if x >= 350]) for j in range(0,len(Is))]

unpicklefile = open('steadystate_0_0_0.sav', 'r')
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()
times_control = unpickledlist[1][5]
Vsoma_control = unpickledlist[2][5]

unpicklefile = open('DCshortthreshs_0_0_0.sav', 'r')
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()
threshIs_control = unpickledlist[0][5]

variants = [[8,0,0],[8,4,0],[12,0,0],[12,2,0]]

close("all")
f,axarr = plt.subplots(2,len(variants))
iters = [0,2,5,6,8]

timesAll = []
VsomaAll = []
spikeFreqsAll = []
threshIsAll = []

for ivar in range(0,len(variants)):
      igene = variants[ivar][0]
      imut = variants[ivar][1]

      nVals = len(MT[igene][imut])*[0]
      thesemutvars = []
      theseCoeffs = theseCoeffsAll[igene][imut]
      for imutvar in range(0,len(MT[igene][imut])):
        thesemutvars.append(MT[igene][imut][imutvar][0])
        if type(MT[igene][imut][imutvar][1]) is int or type(MT[igene][imut][imutvar][1]) is float:
          MT[igene][imut][imutvar][1] = [MT[igene][imut][imutvar][1]]
        nVals[imutvar] = len(MT[igene][imut][imutvar][1])
      cumprodnVals = cumprod(nVals)
      allmutvars = cumprodnVals[len(MT[igene][imut])-1]*[thesemutvars]
      allmutvals = []
      for iallmutval in range(0,cumprodnVals[len(MT[igene][imut])-1]):
        allmutvals.append([0]*len(thesemutvars))
      for iallmutval in range(0,cumprodnVals[len(MT[igene][imut])-1]):
        for imutvar in range(0,len(MT[igene][imut])):
          if imutvar==0:
            allmutvals[iallmutval][imutvar] = MT[igene][imut][imutvar][1][iallmutval%nVals[imutvar]]
          else:
            allmutvals[iallmutval][imutvar] = MT[igene][imut][imutvar][1][(iallmutval/cumprodnVals[imutvar-1])%nVals[imutvar]]
   
      iallmutval = variants[ivar][2]

      unpicklefile = open('ifcurves_'+str(igene)+'_'+str(imut)+'_'+str(iallmutval)+'.sav', 'r')
      unpickledlist = pickle.load(unpicklefile)
      unpicklefile.close()
      ISIsThisMutVal = unpickledlist[0]
      spTimesThisMutVal = unpickledlist[1+ispDef]

      unpicklefile = open('steadystate_'+str(igene)+'_'+str(imut)+'_'+str(iallmutval)+'.sav', 'r')
      unpickledlist = pickle.load(unpicklefile)
      unpicklefile.close()
      timesThisVal = unpickledlist[1]
      VsomaThisVal = unpickledlist[2]

      unpicklefile = open('DCshortthreshs_'+str(igene)+'_'+str(imut)+'_'+str(iallmutval)+'.sav', 'r')
      unpickledlist = pickle.load(unpicklefile)
      unpicklefile.close()
      threshIsThisVal = unpickledlist[0]

      mutText = ""
      for imutvar in range(0,len(MT[igene][imut])):
        if imutvar > 0 and imutvar%2==0:
          mutText = mutText+"\n"
        mutvars = allmutvars[iallmutval][imutvar]
        mutvals = allmutvals[iallmutval][imutvar]
        if type(mutvars) is str:
          mutvars = [mutvars]
        mutText = mutText + str(mutvars) + ": "
        for kmutvar in range(0,len(mutvars)):
          mutvar = mutvars[kmutvar]
          if (mutvar.find('off') > -1 and mutvar.find('offc') < 0) or mutvar.find('eh') > -1:
            if mutvals >= 0 and kmutvar==0:
              mutText = mutText + "+" + str(mutvals) +" mV"
            elif kmutvar==0:
              mutText = mutText  + str(mutvals) +" mV"
          else:
            if kmutvar==0:
              mutText = mutText + "*" + str(mutvals)
          if kmutvar < len(mutvars)-1:
            mutText = mutText + ", "
      print mutText

      nSpsThis = []
      spikeFreqsThisVal = []
      for iiter in range(0,len(iters)):
        iter = iters[iiter]
        if iter == 5:
          continue
        nSpikes = [sum([1 for x in spTimesThisMutVal[iiter][j] if x >= 350]) for j in range(0,len(Is))]
        spikeFreqsThisVal.append([x/15.65 for x in nSpikes])
        axarr[0,ivar].plot(Is, [x/15.65 for x in nSpikes], styles[iter],color=cols[iter])
        axarr[1,ivar].plot(timesThisVal[iiter], VsomaThisVal[iiter], styles[iter],color=cols[iter])

      spikeFreqs_control = [x/15.65 for x in nSpikes_control]
      axarr[0,ivar].plot(Is, [x/15.65 for x in nSpikes_control], styles[iter],color=col_control)
      axarr[1,ivar].plot(times_control, Vsoma_control, styles[iter],color=col_control)
      axarr[0,ivar].set_ylim([0,50])
      axarr[1,ivar].set_xlim([190,600])
      axarr[0,ivar].set_xlabel('I (nA)')
      axarr[1,ivar].set_xlabel('t (ms)')
      if ivar==0:
        axarr[0,ivar].set_ylabel('F (Hz)')
        axarr[0,ivar].set_ylabel('V_m (mV)')
      #timesAll.append(timesThisVal[:])
      #VsomaAll.append(VsomaThisVal[:])
      timesAll.append([timesThisVal[i] for i in [0,1,3,4]])
      VsomaAll.append([VsomaThisVal[i] for i in [0,1,3,4]])
      threshIsAll.append([threshIsThisVal[i] for i in [0,1,3,4]])
      spikeFreqsAll.append(spikeFreqsThisVal[:])

f.suptitle(mutText)
f.savefig("fig1.eps")

scipy.io.savemat('fig1_curves.mat', {'Is': Is, 'spikeFreqs': spikeFreqsAll, 'times': timesAll, 'Vsoma': VsomaAll, 'threshIs': threshIsAll,
                                     'spikeFreqs_control': spikeFreqs_control, 'times_control': times_control, 'Vsoma_control': Vsoma_control, 'threshI_control': threshIs_control,
                                     'variants': variants})

picklelist = [Is, spikeFreqsAll, timesAll, VsomaAll, threshIsAll, spikeFreqs_control, times_control, Vsoma_control, threshIs_control, variants]
file = open('fig1_curves.sav', 'w')
pickle.dump(picklelist,file)
file.close()


