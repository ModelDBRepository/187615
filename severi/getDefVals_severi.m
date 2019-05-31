function [params, varargout] =...
  severi_fantini_charawi_difrancesco_2012_default()
  % Default values for ODE model: severi_fantini_charawi_difrancesco_2012
  % ---------------------------------------------------------------------
  %;
  % params = severi_fantini_charawi_difrancesco_2012_default();
  % [params, ic] = severi_fantini_charawi_difrancesco_2012_default();
  % [params, ic, state_names] =
  % severi_fantini_charawi_difrancesco_2012_default();

  if nargout < 1 || nargout > 3
    error('Expected 1-3 output arguments.');
  end

  % --- Default parameters values --- 

  % --- Voltage clamp ---
  params.V_holding = -45;
  params.V_test = -35;
  params.t_holding = 0.5;
  params.t_test = 0.5;

  % --- Rate modulation experiments ---
  params.ACh = 0;
  params.BAPTA_10_mM = 0;
  params.Cs_5_mM = 0;
  params.Iso_1_uM = 0;
  params.Iva_3_uM = 0;

  % --- I cal fca gate ---
  params.Km_fCa = 0.00035;
  params.alpha_fCa = 0.01;

  % --- Ca sr release ---
  params.EC50_SR = 0.45;
  params.HSR = 2.5;
  params.MaxSR = 15;
  params.MinSR = 1;
  params.kiCa = 500;
  params.kim = 5;
  params.koCa = 10000;
  params.kom = 60;
  params.ks = 250000000;

  % --- Ca intracellular fluxes ---
  params.K_up = 0.0006;
  params.P_up_basal = 12;
  params.tau_dif_Ca = 4e-05;
  params.tau_tr = 0.04;

  % --- Ca buffering ---
  params.CM_tot = 0.045;
  params.CQ_tot = 10;
  params.Mgi = 2.5;
  params.TC_tot = 0.031;
  params.TMC_tot = 0.062;
  params.kb_CM = 542;
  params.kb_CQ = 445;
  params.kb_TC = 446;
  params.kb_TMC = 7.51;
  params.kb_TMM = 751;
  params.kf_CM = 227700;
  params.kf_CQ = 534;
  params.kf_TC = 88800;
  params.kf_TMC = 227700;
  params.kf_TMM = 2277;

  % --- Cell parameters ---
  params.L_cell = 70;
  params.L_sub = 0.02;
  params.R_cell = 4;
  params.V_i_part = 0.46;
  params.V_jsr_part = 0.0012;
  params.V_nsr_part = 0.0116;

  % --- Ionic values ---
  params.Cao = 1.8;
  params.Ki = 140;
  params.Ko = 5.4;
  params.Nao = 140;

  % --- I f ---
  params.Km_f = 45;

  % --- I nak ---
  params.Km_Kp = 1.4;
  params.Km_Nap = 14;
  params.i_NaK_max = 0.063;

  % --- I naca ---
  params.K1ni = 395.3;
  params.K1no = 1628;
  params.K2ni = 2.289;
  params.K2no = 561.4;
  params.K3ni = 26.44;
  params.K3no = 4.663;
  params.K_NaCa = 4;
  params.Kci = 0.0207;
  params.Kcni = 26.44;
  params.Kco = 3.663;
  params.Qci = 0.1369;
  params.Qco = 0;
  params.Qn = 0.4315;

  % --- I na ---
  params.g_Na = 0.0125;

  % --- I na m gate ---
  params.delta_m = 1e-05;

  % --- I cal ---
  params.P_CaL = 0.2;

  % --- I cat ---
  params.P_CaT = 0.02;

  % --- Ca dynamics ---
  params.T_Ca = 6.928;
  params.kbBAPTA = 119.38;
  params.kfBAPTA = 940000;

  % --- I to ---
  params.g_to = 0.002;

  % --- I kr ---
  params.g_Kr = 0.0021637;

  % --- I ks n gate ---
  params.shift = 0;

  % --- I kach ---
  params.g_KACh = 0.00864;

  % --- Membrane ---
  params.C = 3.2e-05;
  params.F = 96485.3415;
  params.R = 8314.472;
  params.T = 310;
  params.clamp_mode = 0;

  
  params.offh_f = -52.5;
  params.sloh_f = 9;
  params.tauh_f = 0.7166529;
  params.offm1_Na = -41.0;
  params.offm2_Na = -66.0;
  params.slom1_Na = 10.0;  
  params.slom2_Na = 1/0.056;
  params.taum_Na = 1;
  params.offh_Na = -75;
  params.sloh1_Na = 8;
  params.sloh2_Na = 10;
  params.tauh_Na = 1;
  params.offm_CaL = -20.3;
  params.slom_CaL = 4.2;
  params.taum_CaL = 0.001;
  params.offh_CaL = -37.4;
  params.sloh_CaL = 5.3;
  params.tauh_CaL = 0.001;
  params.offm_CaT = -38.3;
  params.slom_CaT = 5.5;
  params.taum_CaT = 0.001;
  params.offh_CaT = -58.7;
  params.sloh_CaT = 3.8;
  params.tauh_CaT = 1.0;
  params.offm_Kr = -14.8;
  params.slom_Kr = 8.5;
  params.taum1_Kr = 0.84655354;
  params.taum2_Kr = 1.0;
  params.offh_Kr = -28.6;
  params.sloh_Kr = 17.1;
  params.tauh_Kr = 1.0;
  params.offm1_Ks = 40;
  params.offm2_Ks = 0;
  params.slom1_Ks = 12;
  params.slom2_Ks = 45;
  params.taum_Ks = 1.0;
    
  
  if nargout == 2

    % --- Default initial state values --- 
    x0 = zeros(33, 1);

    % --- I cal fca gate ---
    x0(1) = 0.69799854326;

    % --- Ca sr release ---
    x0(2) = 7.86181717518e-08;
    x0(3) = 1.7340201253e-07;
    x0(4) = 0.211148145513;
    x0(5) = 0.912317231017;

    % --- Ca buffering ---
    x0(6) = 0.0373817991524;
    x0(7) = 0.054381370046;
    x0(8) = 0.299624275429;
    x0(9) = 0.0180519400676;
    x0(10) = 0.281244308217;
    x0(11) = 0.501049376634;

    % --- I f y gate ---
    x0(12) = 0.181334538702;

    % --- I na m gate ---
    x0(13) = 0.440131579216;

    % --- I na h gate ---
    x0(14) = 1.36769401401e-05;

    % --- I cal dl gate ---
    x0(15) = 0;

    % --- I cal fl gate ---
    x0(16) = 0.497133507286;

    % --- I cat dt gate ---
    x0(17) = 0;

    % --- I cat ft gate ---
    x0(18) = 0;

    % --- Ca dynamics ---
    x0(19) = 0.316762674605;
    x0(20) = 1.05386465081;
    x0(21) = 1e-05;
    x0(22) = 1e-05;
    x0(23) = 0;
    x0(24) = 0;

    % --- I to q gate ---
    x0(25) = 0.506139850982;

    % --- I to r gate ---
    x0(26) = 0.0144605370598;

    % --- I kr pa gate ---
    x0(27) = 0.0990510403259;
    x0(28) = 0.322999177803;

    % --- I kr pi gate ---
    x0(29) = 0.705410877259;

    % --- I ks n gate ---
    x0(30) = 0;

    % --- I kach a gate ---
    x0(31) = 0;

    % --- Membrane ---
    x0(32) = -52;

    % --- Nai concentration ---
    x0(33) = 7.5;
    varargout(1) = {x0};
  end

  if nargout == 3

    % --- State names --- 
    state_names = cell(33, 1);

    % --- I cal fca gate ---
    state_names{1} = 'fCa';

    % --- Ca sr release ---
    state_names{2} = 'I';
    state_names{3} = 'O';
    state_names{4} = 'RI';
    state_names{5} = 'R_Ca';

    % --- Ca buffering ---
    state_names{6} = 'fCMi';
    state_names{7} = 'fCMs';
    state_names{8} = 'fCQ';
    state_names{9} = 'fTC';
    state_names{10} = 'fTMC';
    state_names{11} = 'fTMM';

    % --- I f y gate ---
    state_names{12} = 'y';

    % --- I na m gate ---
    state_names{13} = 'm';

    % --- I na h gate ---
    state_names{14} = 'h';

    % --- I cal dl gate ---
    state_names{15} = 'dL';

    % --- I cal fl gate ---
    state_names{16} = 'fL';

    % --- I cat dt gate ---
    state_names{17} = 'dT';

    % --- I cat ft gate ---
    state_names{18} = 'fT';

    % --- Ca dynamics ---
    state_names{19} = 'Ca_jsr';
    state_names{20} = 'Ca_nsr';
    state_names{21} = 'Ca_sub';
    state_names{22} = 'Cai';
    state_names{23} = 'fBAPTA';
    state_names{24} = 'fBAPTA_sub';

    % --- I to q gate ---
    state_names{25} = 'q';

    % --- I to r gate ---
    state_names{26} = 'r';

    % --- I kr pa gate ---
    state_names{27} = 'paF';
    state_names{28} = 'paS';

    % --- I kr pi gate ---
    state_names{29} = 'piy';

    % --- I ks n gate ---
    state_names{30} = 'n';

    % --- I kach a gate ---
    state_names{31} = 'a';

    % --- Membrane ---
    state_names{32} = 'V_ode';

    % --- Nai concentration ---
    state_names{33} = 'Nai_';
    varargout(2) = {state_names};
  end
end