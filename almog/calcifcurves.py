#hoc-code based on implementation https://senselab.med.yale.edu/modeldb/showModel.cshtml?model=151825
#(CC BY 3.0) Tuomo Maki-Marttunen
from neuron import h
import matplotlib
matplotlib.use('Agg')
import numpy
from pylab import *
import mytools
import pickle
import sys


v0 = -62
ca0 = 0.0001
proximalpoint = 400
distalpoint = 620
#distalpoint = 960
BACdt = 5.0
Is = [0.7+0.005*x for x in range(0,41)]
coeffCoeffs = [[0.25,0],[0.125,0],[0.5,0],[0.5,1.0/3],[0.5,2.0/3],[0.5,1.0],[-0.25,0],[-0.125,0],[-0.5,0]]

import mutation_stuff
MT = mutation_stuff.getMT()
defVals = mutation_stuff.getdefvals()
keyList = defVals.keys()
mySuffixes = mutation_stuff.getsuffixes()
mySuffixExceptions = mutation_stuff.getsuffixexceptions()

unpicklefile = open('scalings.sav', 'r')
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()

theseCoeffsAll = unpickledlist[0]
theseMutValsAll = unpickledlist[2]

for icell in range(0,1):

  h("""
load_file("myrun.hoc")
objref cvode
cvode = new CVode()
cvode.active(1)
cvode.atol(0.001)

access a_soma
objref st1,syn1, sl
a_soma st1 = new IClamp(0.5)

double siteVec[2]
sl = new List()
sl=locateSites("apic",620)
maxdiam = 0
for(i=0;i<sl.count();i+=1){
  dd1 = sl.o[i].x[1]
  dd = apic[sl.o[i].x[0]].diam(dd1)
  if (dd > maxdiam) {
    j = i
    maxdiam = dd
  }
}
siteVec[0] = sl.o[j].x[0]
siteVec[1] = sl.o[j].x[1]
apic[siteVec[0]] syn1 = new AlphaSynapse(siteVec[1])
//apic[41] syn1 = new AlphaSynapse(0.5)

syn1.onset = 3400
syn1.tau = 3
syn1.gmax = 0.0
syn1.e = 50

objref vsoma, vdend, tvec
vsoma = new Vector()
vdend = new Vector()
tvec = new Vector()
a_soma cvode.record(&v(0.5),vsoma,tvec)
apic[siteVec[0]] cvode.record(&v(siteVec[1]),vdend,tvec)

v_init = -62
dt = 0.025
tstop = 1000
""")

  styles = ['g-','g-','g-','g-','g-','g-','g-','g-','g-']
  #cols = ['#00aaaa','#11cc44','#55ee00','#bbaa00','#ee6600','#ff0000', '#aa00aa','#772277','#333333']
  cols = ['#666666','#012345','#aa00aa','#bbaa00','#ee6600','#ff0000', '#00aaaa','#772277','#00cc00']
  
  counter = -1
  for igene in range(0,len(MT)):
   for imut in range(0,len(MT[igene])):
    spTimesThisMut = []
    spTimesThisMut2 = []
    ISIs_thismut = []
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
  
    for iallmutval in range(0,cumprodnVals[len(MT[igene][imut])-1]):
      counter = counter + 1
      if len(sys.argv) > 1 and int(float(sys.argv[1])) != counter:
        continue
      mutval = allmutvals[iallmutval]
      nextCoeffs = [0.0,2.0,1.0]
      spTimesThisVal = []
      spTimesThisVal2 = []
      ISIs_thismutval = []
  
      close("all")
      f, axarr = plt.subplots(1, 1)
      axarr.set_position([0.13, 0.1, 0.85, 0.67])
      for iter in [0, 2, 5, 6, 8, -1]:
        if iter >= 0:
          thisCoeff = coeffCoeffs[iter][0]*theseCoeffs[iallmutval] + coeffCoeffs[iter][1]*(1.0 - 0.5*theseCoeffs[iallmutval])
        else:
          thisCoeff = 0
        if iter == -1 and (igene > 0 or imut > 0 or iallmutval > 0):
          continue # do the control only once!
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
            if (mutvars[kmutvar].find('off') > -1 and mutvars[kmutvar].find('offc') < 0) or mutvars[kmutvar].find('eh') > -1:
              newVal =  defVals[mutvar]+mutvals*thisCoeff
              if mutvals >= 0 and kmutvar==0:
                mutText = mutText + "+" + str(mutvals) +" mV"
              elif kmutvar==0:
                mutText = mutText  + str(mutvals) +" mV"
            else:
              newVal = defVals[mutvar]*(mutvals**thisCoeff)
              if kmutvar==0:
                mutText = mutText + "*" + str(mutvals)
            if kmutvar < len(mutvars)-1:
              mutText = mutText + ", "
            mySuffix = mutvars[kmutvar][mutvars[kmutvar].find('_')+1:len(mutvars[kmutvar])]
            mySuffixInd = next((i for i,x in enumerate(mySuffixes) if x.find(mySuffix) > -1))
            isException = 0
            for jsuffe in range(0,len(mySuffixExceptions[mySuffixInd])):
              if mySuffixExceptions[mySuffixInd][jsuffe][0].find(mutvars[kmutvar]) > -1:
                isException = 1
                exceptionInd = jsuffe

            if not isException:
              print ("""forall if(ismembrane(\""""+mySuffix+"""\")) """+mutvars[kmutvar]+""" = """+str(newVal))
              h("""forall if(ismembrane(\""""+mySuffix+"""\")) """+mutvars[kmutvar]+""" = """+str(newVal))
            else:
              print ("""forall if(ismembrane(\""""+mySuffix+"""\")) """+mySuffixExceptions[isuffix][j][1]+""" = """+str(newVal))
              h("""forall if(ismembrane(\""""+mySuffix+"""\")) """+mySuffixExceptions[isuffix][j][1]+""" = """+str(newVal))

        print mutText
        thisCa = h.a_soma.cainf_cad

        spTimesThisCoeff = []
        spTimesThisCoeff2 = []
        ISIs = len(Is)*[0.0]
        nSpikes = []
        for iI in range(0,len(Is)):
          tstop = 16000.0
          squareAmp = Is[iI]
          squareDur = 15800.0
          h("""
tstop = """+str(tstop)+"""
v_init = """+str(v0)+"""
cai0_ca_ion = """+str(thisCa)+"""
st1.amp = """+str(squareAmp)+"""
st1.del = 200
st1.dur = """+str(squareDur)+"""
syn1.gmax = 0
syn1.onset = 200 + """+str(BACdt)+""" 
  """)
          h.init()
          h.run()
  
          times=np.array(h.tvec)
          Vsoma=np.array(h.vsoma)
          Vdend=np.array(h.vdend)
          spikes = mytools.spike_times(times,Vsoma,-50,-50)
          spikes2 = mytools.spike_times(times,Vsoma,-50,inf)
          spTimesThisCoeff.append(spikes[:])
          spTimesThisCoeff2.append(spikes2[:])
          nSpikes1 = len(spikes)
          nSpikes2 = sum([1 for x in spikes if x >= 1000.0])
          nSpikes.append(nSpikes2)
  
          if nSpikes1 > 5:
            spts = spikes[len(spikes)-5:len(spikes)]
            ISIs[iI] = mean([y-x for x,y in zip(spts[0:4],spts[1:5])])
          else:
            ISIs[iI] = 1.0e10
  
        #if iter==0:
        #  axarr.plot(Is, [x/15.0 for x in nSpikes_control])
        axarr.plot(Is, [x/15.0 for x in nSpikes], styles[iter],color=cols[iter])
  
        #Restore default values:
        for imutvar in range(0,len(MT[igene][imut])):
          mutvars = allmutvars[iallmutval][imutvar]
          mutvals = allmutvals[iallmutval][imutvar]
          if type(mutvars) is str:
            mutvars = [mutvars]
          for kmutvar in range(0,len(mutvars)):
            newVal = defVals[mutvars[kmutvar]]
            mySuffix = mutvars[kmutvar][mutvars[kmutvar].find('_')+1:len(mutvars[kmutvar])]
            mySuffixInd = next((i for i,x in enumerate(mySuffixes) if x.find(mySuffix) > -1))
            isException = 0
            for jsuffe in range(0,len(mySuffixExceptions[mySuffixInd])):
              if mySuffixExceptions[mySuffixInd][jsuffe][0].find(mutvars[kmutvar]) > -1:
                isException = 1
                exceptionInd = jsuffe
            if not isException:
              h("""forall if(ismembrane(\""""+mySuffix+"""\")) """+mutvars[kmutvar]+""" = """+str(defVals[mutvars[kmutvar]]))
            else:
              h("""forall if(ismembrane(\""""+mySuffix+"""\")) """+mySuffixExceptions[isuffix][j][1]+""" = """+str(defVals[mutvars[kmutvar]]))

        spTimesThisVal.append(spTimesThisCoeff[:])
        spTimesThisVal2.append(spTimesThisCoeff2[:])
        ISIs_thismutval.append(ISIs[:])
      axarr.set_title('f-I curve')
      xlabel('I (nA)')
      ylabel('F (Hz)')
      f.suptitle(mutText)
      f.savefig("ifcurves_"+str(igene)+"_"+str(imut)+"_"+str(iallmutval)+".eps")

      picklelist = [ISIs_thismutval,spTimesThisVal,spTimesThisVal2,MT]
      file = open('ifcurves_'+str(igene)+'_'+str(imut)+'_'+str(iallmutval)+'.sav', 'w')
      pickle.dump(picklelist,file)
      file.close()
