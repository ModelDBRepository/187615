#hoc-code based on the implementation https://senselab.med.yale.edu/modeldb/showModel.cshtml?model=139653
#(CC BY 3.0) Tuomo Maki-Marttunen
from neuron import h
import matplotlib
matplotlib.use('Agg')
import numpy
from pylab import *
import mytools
import pickle
import sys
from os.path import exists

v0 = -80
ca0 = 0.0001
proximalpoint = 400
distalpoint = 620
#distalpoint = 960
BACdt = 5.0
Is = [0.35+0.05*x for x in range(0,22)]
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

threshIsAllAll = []
threshIsAllAll2 = []
ISIs_allAll = []

for icell in range(0,1):
  theseCoeffsAll = theseCoeffsAllAll[icell]
  threshIsAll = []
  threshIsAll2 = []
  ISIs_all = []
  morphology_file = "morphologies/cell"+str(icell+1)+".asc"
  biophys_file = "models/L5PCbiophys3.hoc"
  template_file = "models/L5PCtemplate.hoc"

  h("""
load_file("stdlib.hoc")
load_file("stdrun.hoc")
objref cvode
cvode = new CVode()
cvode.active(1)
cvode.atol(0.001)
load_file("import3d.hoc")
objref L5PC
load_file(\""""+biophys_file+"""\")
load_file(\""""+template_file+"""\")
L5PC = new L5PCtemplate(\""""+morphology_file+"""\")
objref st1
st1 = new IClamp(0.5)
L5PC.soma st1
L5PC.distribute_channels("apic","gIhbar_Ih",2,-0.8696,3.6161,0.0,1.0*2.0870,0.0002)
L5PC.distribute_channels("apic","gCa_HVAbar_Ca_HVA",3,1.0,0.1,685.0,885.0,1.0*0.000555)
L5PC.distribute_channels("apic","gCa_LVAstbar_Ca_LVAst",3,1.0,0.01,685.0,885.0,1.0*0.0187)
objref sl,st2,ns,syn1,con1,isyn, tvec
isyn = new Vector()
tvec = new Vector()
sl = new List()
double siteVec[2]
sl = L5PC.locateSites("apic","""+str(distalpoint)+""")
maxdiam = 0
for(i=0;i<sl.count();i+=1){
  dd1 = sl.o[i].x[1]
  dd = L5PC.apic[sl.o[i].x[0]].diam(dd1)
  if (dd > maxdiam) {
    j = i
    maxdiam = dd
  }
}
siteVec[0] = sl.o[j].x[0]
siteVec[1] = sl.o[j].x[1]
print "distalpoint gCa_HVA: ", L5PC.apic[siteVec[0]].gCa_HVAbar_Ca_HVA
print "distalpoint gCa_LVA: ", L5PC.apic[siteVec[0]].gCa_LVAstbar_Ca_LVAst
access L5PC.apic[siteVec[0]]
st2 = new IClamp(siteVec[1])
st2.amp = 0
L5PC.apic[siteVec[0]] {
  st2
  syn1 = new epsp(siteVec[1])
  syn1.tau0 = 0.5
  syn1.tau1 = 5
  syn1.onset = 200 + """+str(BACdt)+""" 
  cvode.record(&syn1.i,isyn,tvec)
}
objref vsoma, vdend, recSite, vdend2, isoma, cadend, cadend2, casoma
vsoma = new Vector()
casoma = new Vector()
vdend = new Vector()
cadend = new Vector()
vdend2 = new Vector()
cadend2 = new Vector()
access L5PC.soma
cvode.record(&v(0.5),vsoma,tvec)
cvode.record(&cai(0.5),casoma,tvec)
access L5PC.apic[siteVec[0]]
cvode.record(&v(siteVec[1]),vdend,tvec)
cvode.record(&cai(siteVec[1]),cadend,tvec)
sl = new List()
sl = L5PC.locateSites("apic","""+str(proximalpoint)+""")
maxdiam = 0
for(i=0;i<sl.count();i+=1){
  dd1 = sl.o[i].x[1]
  dd = L5PC.apic[sl.o[i].x[0]].diam(dd1)
  if (dd > maxdiam) {
    j = i
    maxdiam = dd
  }
}
siteVec[0] = sl.o[j].x[0]
siteVec[1] = sl.o[j].x[1]
access L5PC.apic[siteVec[0]]
recSite = new IClamp(siteVec[1])
recSite.amp = 0
L5PC.apic[siteVec[0]] {
        recSite
}
access L5PC.apic[siteVec[0]]
cvode.record(&v(siteVec[1]),vdend2,tvec)
cvode.record(&cai(siteVec[1]),cadend2,tvec)
access L5PC.soma
isoma = new Vector()
cvode.record(&st1.i,isoma,tvec)
""")

  styles = ['g-','g-','g-','g-','g-','g-','g-','g-','g-']
  #cols = ['#00aaaa','#11cc44','#55ee00','#bbaa00','#ee6600','#ff0000', '#aa00aa','#772277','#333333']
  cols = ['#666666','#012345','#aa00aa','#bbaa00','#ee6600','#ff0000', '#00aaaa','#772277','#00cc00']
  
  counter = -1
  for igene in range(0,len(MT)):
   threshIsThisGene = []
   threshIsThisGene2 = []
   ISIs_thisgene = []
   for imut in range(0,len(MT[igene])):
    threshIsThisMut = []
    threshIsThisMut2 = []
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
      threshIsThisVal = []
      threshIsThisVal2 = []
      ISIs_thismutval = []
      if exists('DCshortthreshs_cs'+str(icell)+'_'+str(igene)+'_'+str(imut)+'_'+str(iallmutval)+'.sav'):
        continue
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
            if mutvar.find('offm') > -1 or mutvar.find('offh') > -1 or mutvar.find('ehcn') > -1:
              newVal =  [x+mutvals*thisCoeff for x in defVals[mutvar]]
              if mutvals >= 0 and kmutvar==0:
                mutText = mutText + "+" + str(mutvals) +" mV"
              elif kmutvar==0:
                mutText = mutText  + str(mutvals) +" mV"
            else:
              newVal =  [x*(mutvals**thisCoeff) for x in defVals[mutvar]]
              if kmutvar==0:
                mutText = mutText + "*" + str(mutvals)
            if kmutvar < len(mutvars)-1:
              mutText = mutText + ", "
            #else:                                                                                                                                                                       
            #  mutText = mutText + "\n"                                                                                                                                                  
            if mutvar.find('_Ih') > -1:
              updateThese = [1,1,1]
            elif mutvar.find('_Ca_HVA') > -1 or mutvar.find('_Ca_LVAst') > -1 or mutvar.find('_SKv3.1') > -1 or mutvar.find('_Ca_HVA') > -1 or mutvar.find('_SK_E2') > -1 or mutvar.find\
  ('_NaTa_t') > -1 or mutvar.find('_CaDynamics_E2') > -1:
              updateThese = [1,1,0]
            elif mutvar.find('_K_Pst') > -1 or mutvar.find('_K_Tst') > -1 or mutvar.find('_Nap_Et2') > -1:
              updateThese = [1,0,0]
            elif mutvar.find('_Im') > -1:
              updateThese = [0,1,0]
            else:
              print "Error: str=" + str(mutvar)
              updatedThese = [0,0,0]
            for iupdated in range(0,3):
              if updateThese[iupdated]:
                print """forsec L5PC."""+str(updatedVars[iupdated])+""" {                                                                                                                
  """+mutvar+""" = """+str(newVal[whichDefVal[iupdated]])+"""                                                                                                                            
  }"""
                h("""forsec L5PC."""+str(updatedVars[iupdated])+""" {                                                                                                                    
  """+mutvar+""" = """+str(newVal[whichDefVal[iupdated]])+"""                                                                                                                            
  }""")
        print mutText
        ISIs = len(Is)*[0.0]
        nSpikes = []

        ITERS = 30
        nextIs = [0.0,0.5,0.25]

        for iITER in range(0,ITERS+2):
          tstop = 4000.0
          squareAmp = nextIs[min(iITER,2)]
          squareDur = 3800.0
          epsp_Imax = 0.0
          h("""
  tstop = """+str(tstop)+"""
  v_init = """+str(v0)+"""
  cai0_ca_ion = """+str(ca0)+"""
  st1.amp = """+str(squareAmp)+"""
  st1.del = 200
  st1.dur = """+str(squareDur)+"""
  syn1.imax = """+str(epsp_Imax)+"""
  syn1.onset = 200 + """+str(BACdt)+""" 
  """)
          h.init()
          h.run()
  
          times=np.array(h.tvec)
          Vsoma=np.array(h.vsoma)
          Vdend=np.array(h.vdend)
          spikes = mytools.spike_times(times,Vsoma,-35,-45)
          nSpikes2 = sum([1 for x in spikes if x >= 0.0])
  
          isChanged = nSpikes2 > 0
          print str(nSpikes2)+" spikes, I="+str(squareAmp)
          if iITER==0 and isChanged:
            print "Even zero amplitude causes spiking!! igene="+str(igene)+", imut="+str(imut)+", iallmutval="+str(iallmutval)
            continue
          if iITER==1 and not isChanged:
            print "Even large current does not cause spiking!! igene="+str(igene)+", imut="+str(imut)+", iallmutval="+str(iallmutval)
            continue
          if iITER>=2 and iITER < ITERS+2:
            if isChanged:
              nextIs = [nextIs[0],nextIs[2],0.5*nextIs[0]+0.5*nextIs[2]]
            else:
              nextIs = [nextIs[2],nextIs[1],0.5*nextIs[1]+0.5*nextIs[2]]

  
        #Print the parameters and their default values:
        for idefval in range(0,len(defVals.keys())):
          thisdefval = defVals.keys()[idefval]
          if thisdefval.find('_Im') > -1:
            h('print "L5PC.apic[0].'+thisdefval+' = ", L5PC.apic[0].'+thisdefval+', "Default = ", '+str(defVals[thisdefval][1]))
            #) #+" (def="+str(defVals[thisdefval])+")"
          else:
            h('print "L5PC.soma[0].'+thisdefval+' = ", L5PC.soma[0].'+thisdefval+', "Default = ", '+str(defVals[thisdefval][0]))
            #h('print L5PC.soma[0]."+thisdefval) #+" (def="+str(defVals[thisdefval])+")"      
  
        #Restore default values:
        for imutvar in range(0,len(MT[igene][imut])):
          mutvars = allmutvars[iallmutval][imutvar]
          mutvals = allmutvals[iallmutval][imutvar]
          if type(mutvars) is str:
            mutvars = [mutvars]
          for kmutvar in range(0,len(mutvars)):
            mutvar = mutvars[kmutvar]
            newVal = defVals[mutvar]
            if mutvar.find('_Ih') > -1:
              updateThese = [1,1,1]
            elif mutvar.find('_Ca_HVA') > -1 or mutvar.find('_Ca_LVAst') > -1 or mutvar.find('_SKv3.1') > -1 or mutvar.find('_Ca_HVA') > -1 or mutvar.find('_SK_E2') > -1 or mutvar.find('_NaTa_t') > -1 or mutvar.find('_CaDynamics_E2') > -1:
              updateThese = [1,1,0]
            elif mutvar.find('_K_Pst') > -1 or mutvar.find('_K_Tst') > -1 or mutvar.find('_Nap_Et2') > -1:
              updateThese = [1,0,0]
            elif mutvar.find('_Im') > -1:
              updateThese = [0,1,0]
            else:
              print "Error: str=" + str(mutvar)
              updatedThese = [0,0,0]
            for iupdated in range(0,3):
              if updateThese[iupdated]:
                print """forsec L5PC."""+str(updatedVars[iupdated])+""" {
  """+mutvar+""" = """+str(newVal[whichDefVal[iupdated]])+"""
  }"""
                h("""forsec L5PC."""+str(updatedVars[iupdated])+""" {
  """+mutvar+""" = """+str(newVal[whichDefVal[iupdated]])+"""
  }""")
        threshIsThisVal.append(nextIs[2])
      xlabel('I (nA)')
      ylabel('F (Hz)')

      threshIsThisMut.append(threshIsThisVal[:])
      picklelist = [threshIsThisVal,MT]
      file = open('DCshortthreshs_cs'+str(icell)+'_'+str(igene)+'_'+str(imut)+'_'+str(iallmutval)+'.sav', 'w')
      pickle.dump(picklelist,file)
      file.close()
    threshIsThisGene.append(threshIsThisMut[:])
   threshIsAll.append(threshIsThisGene[:])
  
  #picklelist = [ISIs_all,threshIsAll,threshIsAll2,MT]
  #file = open('ifcurves_mut.sav', 'w')
  #pickle.dump(picklelist,file)
  #file.close()
  threshIsAllAll.append(threshIsAll)

#picklelist = [ISIs_allAll,threshIsAllAll,threshIsAllAll2,MT]
#file = open('DCthreshs.sav', 'w')
#pickle.dump(picklelist,file)
#file.close()
  
