function A = getDefVals_kharche()
  A = struct;
  A.offm_st = -67.0;
  A.offm_CaT = -26.0;
  A.slom_CaT = 6.0;
  A.taum_CaT = 1.0;
  A.offh_CaT = -61.7;
  A.sloh_CaT = 5.6;
  A.tauh_CaT = 1.0;
  A.offm_Kr = -21.173694;
  A.slom_Kr = 9.757086;
  A.taum_Kr = 0.699821;
  A.offh_Kr = -16.758474;
  A.sloh_Kr = 19.0;
  A.tauh1_Kr = 0.2;
  A.tauh2_Kr = 0.9;
  A.offm_Ks = 20.876040;
  A.slom_Ks = 11.852723;
  A.taum_Ks = 1000.0;
  A.offm_CaL12 = -3.0;
  A.slom_CaL12 = 5.0;
  A.offh_CaL12 = -36;
  A.sloh_CaL12 = 4.6;
  A.offm_CaL13 = -13.5;
  A.slom_CaL13 = 6.0;
  A.offh_CaL13 = -35;
  A.sloh_CaL13 = 7.3;
  A.taum_CaL = 2000;
  A.tauh_CaL = 1.0;
  A.offm_NaTTXS = -31.097331;
  A.slom_NaTTXS = 5.0;
  A.offh_NaTTXS = -56.0;
  A.sloh_NaTTXS = 3.0;
  A.taum_NaTTXS = 1000.0;
  A.tauh_NaTTXS = 1000.0;
  A.tauj_NaTTXS = 1000.0;
  A.offm_NaTTXR = -45.213705;
  A.slom_NaTTXR = 7.219547;
  A.offh_NaTTXR = -62.578120;
  A.sloh_NaTTXR = 6.084036;
  A.taum_NaTTXR = 1000.0;
  A.tauh_NaTTXR = 1000.0;
  A.tauj_NaTTXR = 1000.0;
  A.offh_f = -106.8;
  A.sloh_f = 16.3;
  A.tauh_f = 1.5049;
  A.Pup_SERCA = 0.04;

  A.g_st = 0.00006;
  A.g_bNa = 0.0001215;
  A.g_bCa = 0.000015;
  A.g_bK = 0.0000025;
  A.g_K1 = 0.229*0.0039228*0.9;
  A.g_Ks = 0.000299;
  A.g_sus = 0.00039060;
  A.g_NaTTXS = 0.1*5.925e-05;
  A.g_NaTTXR = 0.1*5.925e-05;
  A.g_CaL12 = 0.0010*4.0*1.5;
  A.g_CaL13 = 0.0030*4.0*1.5;
  A.g_CaT = 0.75*0.01862; 
  A.g_f = 0.0057; 
  A.g_Kr = 0.8*0.002955;
  A.g_to = 0.00492;
  A.i_NaK = 0.077;
  A.k_NaCa = 5.5;

  A.nao = 140.0;
  A.cao = 1.8;
  A.ko = 5.4;
  
  A.T = 310.5;
  
end