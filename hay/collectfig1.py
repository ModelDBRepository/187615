import matplotlib
matplotlib.use('Agg')
import numpy
from pylab import *
import mytools
import pickle
import sys
import scipy.io

Is = [0.2+0.05*x for x in range(0,25)]
coeffCoeffs = [[0.25,0],[0.125,0],[0.5,0],[0.5,1.0/3],[0.5,2.0/3],[0.5,1.0],[-0.25,0],[-0.125,0],[-0.5,0]]

import mutation_stuff
MT = mutation_stuff.getMT()
defVals = mutation_stuff.getdefvals()
keyList = defVals.keys()
for idefval in range(0,len(keyList)):
  if type(defVals[keyList[idefval]]) is not list:
    defVals[keyList[idefval]] = [defVals[keyList[idefval]], defVals[keyList[idefval]]] #make the dictionary values [somatic, apical]
updatedVars = ['somatic','apical','basal'] # the possible classes of segments that defVals may apply to
whichDefVal = [0,1,0]                      # use the defVal[0] for somatic and basal segments and defVal[1] for apical segments
unpicklefile = open('scalings_cs.sav', 'r')
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()

theseCoeffsAllAll = unpickledlist[0]
theseMutValsAll = unpickledlist[2]

ispDef = 0 # Consider a local maximum above -35mV a spike only if after last spike the membrane potential came below -45mV
#ispDef = 1 # Consider every local maximum above -35mV a spike 

variants = [[8,0,0],[8,4,0],[12,0,0],[12,2,0]]

for icell in range(0,1):
  theseCoeffsAll = theseCoeffsAllAll[icell]

  unpicklefile = open('ifcurves_cs'+str(icell)+'_0_0_0.sav', 'r')
  unpickledlist = pickle.load(unpicklefile)
  unpicklefile.close()
  spTimesThisMutVal = unpickledlist[1+ispDef]
  nSpikes_control = [sum([1 for x in spTimesThisMutVal[5][j] if x >= 500]) for j in range(0,len(Is))]

  unpicklefile = open('steadystate_cs'+str(icell)+'_0_0_0.sav', 'r')
  unpickledlist = pickle.load(unpicklefile)
  unpicklefile.close()
  times_control = unpickledlist[1][5]
  Vsoma_control = unpickledlist[2][5]
  
  unpicklefile = open('DCshortthreshs_cs'+str(icell)+'_0_0_0.sav', 'r')
  unpickledlist = pickle.load(unpicklefile)
  unpicklefile.close()
  threshIs_control = unpickledlist[0][5]

  styles = ['g-','g-','g-','g-','g-','g-','g-','g-','g-']
  #cols = ['#00aaaa','#11cc44','#55ee00','#bbaa00','#ee6600','#ff0000', '#aa00aa','#772277','#333333']
  cols = ['#666666','#012345','#aa00aa','#bbaa00','#ee6600','#ff0000', '#00aaaa','#772277','#00cc00']
  col_control = '#2222ff'  

  counter = -1

  close("all")                               
  f, axarr = plt.subplots(2, len(variants))              
  #axarr.set_position([0.13, 0.1, 0.85, 0.67])

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

  
      unpicklefile = open('ifcurves_cs'+str(icell)+'_'+str(igene)+'_'+str(imut)+'_'+str(iallmutval)+'.sav', 'r')
      unpickledlist = pickle.load(unpicklefile)
      unpicklefile.close()
      ISIsThisMutVal = unpickledlist[0]
      spTimesThisMutVal = unpickledlist[1+ispDef]

      try:
        unpicklefile = open('steadystate_cs'+str(icell)+'_'+str(igene)+'_'+str(imut)+'_'+str(iallmutval)+'.sav', 'r')
      except:
        print 'steadystate_cs'+str(icell)+'_'+str(igene)+'_'+str(imut)+'_'+str(iallmutval)+'.sav not found!'
        unpicklefile = open('steadystate_cs'+str(icell)+'_0_0_0.sav', 'r')
      unpickledlist = pickle.load(unpicklefile)
      unpicklefile.close()
      timesThisVal = unpickledlist[1]
      VsomaThisVal = unpickledlist[2]

      unpicklefile = open('DCshortthreshs_cs'+str(icell)+'_'+str(igene)+'_'+str(imut)+'_'+str(iallmutval)+'.sav', 'r')
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
          if mutvar.find('offm') > -1 or mutvar.find('offh') > -1 or mutvar.find('ehcn') > -1:
            if mutvals >= 0 and kmutvar==0:
              mutText = mutText + "+" + str(mutvals) +" mV"
            elif kmutvar==0:
              mutText = mutText  + str(mutvals) +" mV"
          else:
            if kmutvar==0:
              mutText = mutText + "*" + str(mutvals)
          if kmutvar < len(mutvars)-1:
            mutText = mutText + ", "
      mutText = mutText + " coeff = "+str(theseCoeffs[iallmutval])


      iters = [0, 2, 5, 6, 8]
      spikeFreqsThisVal = []
      for iiter in range(0,len(iters)):
        iter = iters[iiter]
        if iter==5:
          continue
        nSpikes = [sum([1 for x in spTimesThisMutVal[iiter][j] if x >= 500]) for j in range(0,len(Is))]
        spikeFreqsThisVal.append([x/7.5 for x in nSpikes])
        axarr[0,ivar].plot(Is, [x/7.5 for x in nSpikes], styles[iter],color=cols[iter])
        axarr[1,ivar].plot(timesThisVal[iiter], VsomaThisVal[iiter], styles[iter],color=cols[iter])
      spikeFreqs_control = [x/7.5 for x in nSpikes_control]
      axarr[0,ivar].plot(Is, [x/7.5 for x in nSpikes_control], styles[iter],color=col_control)
      axarr[1,ivar].plot(times_control, Vsoma_control, styles[iter],color=col_control)
      #axarr.set_title('I-F curve')
      axarr[0,ivar].set_ylim([0,25])
      axarr[1,ivar].set_xlim([190,600])

      #axarr.set_title('I-F curve')
      axarr[0,ivar].set_ylim([0,25])
      
      axarr[0,ivar].set_xlabel('I (nA)')
      axarr[1,ivar].set_xlabel('t (ms)')
      if ivar==0:
        axarr[0,ivar].set_ylabel('F (Hz)')
        axarr[0,ivar].set_ylabel('V_m (mV)')
      timesAll.append([timesThisVal[i] for i in [0,1,3,4]])
      VsomaAll.append([VsomaThisVal[i] for i in [0,1,3,4]])
      spikeFreqsAll.append(spikeFreqsThisVal[:])
      threshIsAll.append([threshIsThisVal[i] for i in [0,1,3,4]])

  f.suptitle(mutText)
  f.savefig("fig1_cs"+str(icell)+".eps")

  scipy.io.savemat('fig1_curves.mat', {'Is': Is, 'spikeFreqs': spikeFreqsAll, 'times': timesAll, 'Vsoma': VsomaAll, 'threshIs': threshIsAll, 
                                       'spikeFreqs_control': spikeFreqs_control, 'times_control': times_control, 'Vsoma_control': Vsoma_control, 'threshI_control': threshIs_control,
                                       'variants': variants})
  
  picklelist = [Is, spikeFreqsAll, timesAll, VsomaAll, threshIsAll, spikeFreqs_control, times_control, Vsoma_control, threshIs_control, variants]
  file = open('fig1_curves.sav', 'w')
  pickle.dump(picklelist,file)
  file.close()

