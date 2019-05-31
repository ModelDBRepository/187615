#hoc-code based on the implementation https://senselab.med.yale.edu/modeldb/showModel.cshtml?model=139653
#(CC BY 3.0) Tuomo Maki-Marttunen
from neuron import h
import matplotlib
matplotlib.use('Agg')
import numpy
from pylab import *
import mytools
import pickle
import time
import sys
import scipy.io
from os.path import exists

morphology_file = "morphologies/cell1.asc"
biophys_file = "models/L5PCbiophys3.hoc"
template_file = "models/L5PCtemplate.hoc"
v0 = -80
ca0 = 0.0001
proximalpoint = 400
distalpoint = 620
#distalpoint = 960
BACdt = 5.0
fs = 8
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

for icell in range(0,1):
  theseCoeffsAll = theseCoeffsAllAll[icell]
  morphology_file = "morphologies/cell"+str(icell+1)+".asc"
  biophys_file = "models/L5PCbiophys3.hoc"
  template_file = "models/L5PCtemplate.hoc"

  h("""
load_file("stdlib.hoc")
load_file("stdrun.hoc")
objref cvode
cvode = new CVode()
cvode.active(1)
load_file("import3d.hoc")
objref L5PC
load_file(\""""+biophys_file+"""\")
load_file(\""""+template_file+"""\")
L5PC = new L5PCtemplate(\""""+morphology_file+"""\")
objref st1
st1 = new IClamp(0.5)
L5PC.soma st1
forsec L5PC.somatic {
}
forsec L5PC.apical {
}
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
  syn1.onset = 145 + """+str(BACdt)+""" 
  cvode.record(&syn1.i,isyn,tvec)
}
objref vsoma, vdend, recSite, vdend2, isoma, cadend, cadend2, casoma, inatsoma, inapsoma, icahvasoma, icalvasoma, ihsoma, isksoma, inatdend, icahvadend, icalvadend, ihdend, iskdend
objref ikv31soma, iktsoma, ikpsoma, imdend
vsoma = new Vector()
casoma = new Vector()
vdend = new Vector()
cadend = new Vector()
vdend2 = new Vector()
cadend2 = new Vector()
inatsoma = new Vector()
inapsoma = new Vector()
icahvasoma = new Vector()
icalvasoma = new Vector()
ihsoma = new Vector()
isksoma = new Vector()
inatdend = new Vector()
icahvadend = new Vector()
icalvadend = new Vector()
ihdend = new Vector()
iskdend = new Vector()
ikv31soma = new Vector()
iktsoma = new Vector()
ikpsoma = new Vector()
imdend = new Vector()
access L5PC.soma
cvode.record(&v(0.5),vsoma,tvec)
cvode.record(&cai(0.5),casoma,tvec)
cvode.record(&ina_NaTa_t(0.5),inatsoma,tvec)
cvode.record(&ina_Nap_Et2(0.5),inapsoma,tvec)
cvode.record(&ica_Ca_HVA(0.5),icahvasoma,tvec)
cvode.record(&ica_Ca_LVAst(0.5),icalvasoma,tvec)
cvode.record(&ihcn_Ih(0.5),ihsoma,tvec)
cvode.record(&ik_SK_E2(0.5),isksoma,tvec)
cvode.record(&ik_SKv3_1(0.5),ikv31soma,tvec)
cvode.record(&ik_K_Tst(0.5),iktsoma,tvec)
cvode.record(&ik_K_Pst(0.5),ikpsoma,tvec)
access L5PC.apic[siteVec[0]]
cvode.record(&v(siteVec[1]),vdend,tvec)
cvode.record(&cai(siteVec[1]),cadend,tvec)
cvode.record(&ina_NaTa_t(siteVec[1]),inatdend,tvec)
cvode.record(&ica_Ca_HVA(siteVec[1]),icahvadend,tvec)
cvode.record(&ica_Ca_LVAst(siteVec[1]),icalvadend,tvec)
cvode.record(&ihcn_Ih(siteVec[1]),ihdend,tvec)
cvode.record(&ik_SK_E2(siteVec[1]),iskdend,tvec)
cvode.record(&ik_Im(siteVec[1]),imdend,tvec)
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
print "proximalpoint gCa_HVA: ", L5PC.apic[siteVec[0]].gCa_HVAbar_Ca_HVA
print "proximalpoint gCa_LVA: ", L5PC.apic[siteVec[0]].gCa_LVAstbar_Ca_LVAst
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

  counter = -1
  for igene in range(0,len(MT)):
   for imut in range(0,len(MT[igene])):
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
      spTimesThisMutVal = []
      timesThisMutVal = []
      VsomasThisMutVal = []
      CasomasThisMutVal = []
      VdendsThisMutVal = []
      CadendsThisMutVal = []
      InatsomasThisMutVal = []
      InapsomasThisMutVal = []
      IcahvasomasThisMutVal = []
      IcalvasomasThisMutVal = []
      IhsomasThisMutVal = []
      IsksomasThisMutVal = []
      InatdendsThisMutVal = []
      IcahvadendsThisMutVal = []
      IcalvadendsThisMutVal = []
      IhdendsThisMutVal = []
      IskdendsThisMutVal = []
      Ikv31somasThisMutVal = []
      IktsomasThisMutVal = []
      IkpsomasThisMutVal = []
      ImdendsThisMutVal = []
      if exists('steadystate_cs'+str(icell)+'_'+str(igene)+'_'+str(imut)+'_'+str(iallmutval)+'.sav'):
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
              newVal = [x*(mutvals**thisCoeff) for x in defVals[mutvar]]
              if kmutvar==0:
                mutText = mutText + "*" + str(mutvals)
            if kmutvar < len(mutvars)-1:
              mutText = mutText + ", "
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
        print mutText
        tstop = 4000.0
        squareAmp = 1.0
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
  """)
        h.init()
        h.run()
  
        times=np.array(h.tvec)
        Vsoma=np.array(h.vsoma)
        Vdend=np.array(h.vdend)
        Casoma=np.array(h.casoma)
        Cadend=np.array(h.cadend)
        Inatsoma=np.array(h.inatsoma)
        Inapsoma=np.array(h.inapsoma)
        Icahvasoma=np.array(h.icahvasoma)
        Icalvasoma=np.array(h.icalvasoma)
        Ihsoma=np.array(h.ihsoma)
        Isksoma=np.array(h.isksoma)
        Inatdend=np.array(h.inatdend)
        Icahvadend=np.array(h.icahvadend)
        Icalvadend=np.array(h.icalvadend)
        Ihdend=np.array(h.ihdend)
        Iskdend=np.array(h.iskdend)
        Ikv31soma=np.array(h.ikv31soma)
        Iktsoma=np.array(h.iktsoma)
        Ikpsoma=np.array(h.ikpsoma)
        Imdend=np.array(h.imdend)

        spikes = mytools.spike_times(times,Vsoma,-35,-45)
        spTimesThisCoeff = spikes[:]
        nSpikes1 = len(spikes)
  
        if nSpikes1 > 5:
          spts = spikes[len(spikes)-3:len(spikes)]
          istart = next((i for i,x in enumerate(times) if x > spts[0]))
          iend = next((i for i,x in enumerate(times) if x > spts[1]))+4
          nsteps = iend-istart-1
          tdiff = [y-x for x,y in zip(times[istart:iend-1],times[istart+1:iend])]
          cadiff = [y-x for x,y in zip(Casoma[istart:iend-1],Casoma[istart+1:iend])]
          caddiff = [y-x for x,y in zip(Cadend[istart:iend-1],Cadend[istart+1:iend])]
          caderiv1 = [y/x for x,y in zip(tdiff[0:nsteps-1],cadiff[0:nsteps-1])]
          caderiv2 = [y/x for x,y in zip(tdiff[1:nsteps],cadiff[1:nsteps])]
          caderiv = [(x+y)/2.0 for x,y in zip(caderiv1,caderiv2)]
          cadderiv = [y/x for x,y in zip(tdiff,caddiff)]
  
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
        spTimesThisMutVal.append(spTimesThisCoeff[:])
        timesThisMutVal.append(times[:])
        VsomasThisMutVal.append(Vsoma[:])
        CasomasThisMutVal.append(Casoma[:])
        VdendsThisMutVal.append(Vdend[:])
        CadendsThisMutVal.append(Cadend[:])
        InatsomasThisMutVal.append(Inatsoma[:])
        InapsomasThisMutVal.append(Inapsoma[:])
        IcahvasomasThisMutVal.append(Icahvasoma[:])
        IcalvasomasThisMutVal.append(Icalvasoma[:])
        IhsomasThisMutVal.append(Ihsoma[:])
        IsksomasThisMutVal.append(Isksoma[:])
        InatdendsThisMutVal.append(Inatdend[:])
        IcahvadendsThisMutVal.append(Icahvadend[:])
        IcalvadendsThisMutVal.append(Icalvadend[:])
        IhdendsThisMutVal.append(Ihdend[:])
        IskdendsThisMutVal.append(Iskdend[:])
        Ikv31somasThisMutVal.append(Ikv31soma[:])
        IktsomasThisMutVal.append(Iktsoma[:])
        IkpsomasThisMutVal.append(Ikpsoma[:])
        ImdendsThisMutVal.append(Imdend[:])
        if iter==-1:
          picklelist = [theseCoeffsAllAll,times,Vsoma,Casoma,Vdend,Cadend,Inatsoma,Inapsoma,Icahvasoma,Icalvasoma,Ihsoma,Isksoma,Inatdend,Icahvadend,Icalvadend,Ihdend,Iskdend,Ikv31soma,Iktsoma,Ikpsoma,Imdend,MT]
          file = open('steadystate_cs'+str(icell)+'_control.sav', 'w')
          pickle.dump(picklelist,file)
          file.close()
          scipy.io.savemat('steadystate_cs'+str(icell)+'_control.mat', {'times': times, 'Vsoma': Vsoma})
  
      picklelist = [theseCoeffsAllAll,timesThisMutVal,VsomasThisMutVal,CasomasThisMutVal,VdendsThisMutVal,CadendsThisMutVal,
                    InatsomasThisMutVal,InapsomasThisMutVal,IcahvasomasThisMutVal,IcalvasomasThisMutVal,IhsomasThisMutVal,IsksomasThisMutVal,
                    InatdendsThisMutVal,IcahvadendsThisMutVal,IcalvadendsThisMutVal,IhdendsThisMutVal,IskdendsThisMutVal,Ikv31somasThisMutVal,IktsomasThisMutVal,IkpsomasThisMutVal,ImdendsThisMutVal,MT]
      file = open('steadystate_cs'+str(icell)+'_'+str(igene)+'_'+str(imut)+'_'+str(iallmutval)+'.sav', 'w')
      pickle.dump(picklelist,file)
      file.close()

  
