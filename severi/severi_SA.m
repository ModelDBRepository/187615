function [vs,is,ts,peak_times] = severi_SA(T,parChange,PEAK_THR,tolerance)
  if nargin < 1 || isempty(T), T = 1000; end
  if nargin < 2 || isempty(parChange), parChange = struct; end
  if nargin < 3 || isempty(PEAK_THR), PEAK_THR = 0; end %mV; recognise only local maxima that are above PEAK_THR  
  if nargin < 4 || isempty(tolerance), tolerance = [1e-6, 1e-6]; elseif length(tolerance)==1, tolerance = tolerance*[1 1]; end
  
  [p, x] = getDefVals_severi();
  
  if isfield(parChange,'offh_f'), p.offh_f = parChange.offh_f; end
  if isfield(parChange,'sloh_f'), p.sloh_f = parChange.sloh_f; end
  if isfield(parChange,'tauh_f'), p.tauh_f = parChange.tauh_f; end
  if isfield(parChange,'offm1_Na'), p.offm1_Na = parChange.offm1_Na; end
  if isfield(parChange,'offm2_Na'), p.offm2_Na = parChange.offm2_Na; end
  if isfield(parChange,'slom1_Na'), p.slom1_Na = parChange.slom1_Na; end
  if isfield(parChange,'slom2_Na'), p.slom2_Na = parChange.slom2_Na; end
  if isfield(parChange,'taum_Na'), p.taum_Na = parChange.taum_Na; end
  if isfield(parChange,'offh_Na'), p.offh_Na = parChange.offh_Na; end
  if isfield(parChange,'sloh1_Na'), p.sloh1_Na = parChange.sloh1_Na; end
  if isfield(parChange,'sloh2_Na'), p.sloh2_Na = parChange.sloh2_Na; end
  if isfield(parChange,'tauh_Na'), p.tauh_Na = parChange.tauh_Na; end
  if isfield(parChange,'offm_CaL'), p.offm_CaL = parChange.offm_CaL; end
  if isfield(parChange,'slom_CaL'), p.slom_CaL = parChange.slom_CaL; end
  if isfield(parChange,'taum_CaL'), p.taum_CaL = parChange.taum_CaL; end
  if isfield(parChange,'offh_CaL'), p.offh_CaL = parChange.offh_CaL; end
  if isfield(parChange,'sloh_CaL'), p.sloh_CaL = parChange.sloh_CaL; end
  if isfield(parChange,'tauh_CaL'), p.tauh_CaL = parChange.tauh_CaL; end
  if isfield(parChange,'offm_CaT'), p.offm_CaT = parChange.offm_CaT; end
  if isfield(parChange,'slom_CaT'), p.slom_CaT = parChange.slom_CaT; end
  if isfield(parChange,'taum_CaT'), p.taum_CaT = parChange.taum_CaT; end
  if isfield(parChange,'offh_CaT'), p.offh_CaT = parChange.offh_CaT; end
  if isfield(parChange,'sloh_CaT'), p.sloh_CaT = parChange.sloh_CaT; end
  if isfield(parChange,'tauh_CaT'), p.tauh_CaT = parChange.tauh_CaT; end
  if isfield(parChange,'offm_Kr'), p.offm_Kr = parChange.offm_Kr; end
  if isfield(parChange,'slom_Kr'), p.slom_Kr = parChange.slom_Kr; end
  if isfield(parChange,'taum1_Kr'), p.taum1_Kr = parChange.taum1_Kr; end
  if isfield(parChange,'taum2_Kr'), p.taum2_Kr = parChange.taum2_Kr; end
  if isfield(parChange,'offh_Kr'), p.offh_Kr = parChange.offh_Kr; end
  if isfield(parChange,'sloh_Kr'), p.sloh_Kr = parChange.sloh_Kr; end
  if isfield(parChange,'tauh_Kr'), p.tauh_Kr = parChange.tauh_Kr; end
  if isfield(parChange,'offm1_Ks'), p.offm1_Ks = parChange.offm1_Ks; end
  if isfield(parChange,'offm2_Ks'), p.offm2_Ks = parChange.offm2_Ks; end
  if isfield(parChange,'slom1_Ks'), p.slom1_Ks = parChange.slom1_Ks; end
  if isfield(parChange,'slom2_Ks'), p.slom2_Ks = parChange.slom2_Ks; end
  if isfield(parChange,'taum_Ks'), p.taum_Ks = parChange.taum_Ks; end
  if isfield(parChange,'P_up_basal'), p.P_up_basal = parChange.P_up_basal; end

  if isfield(parChange,'Nao'), p.Nao = parChange.Nao; end
  if isfield(parChange,'Ko'), p.Ko = parChange.Ko; end
  if isfield(parChange,'Cao'), p.Cao = parChange.Cao; end
  if isfield(parChange,'Ki'), p.Ki = parChange.Ki; end

  options = odeset('AbsTol',tolerance(1),'RelTol',tolerance(2));
  [ts, xs] = ode15s(@severi_fantini_charawi_difrancesco_2012, [0, T], x, options, p);
  vs = xs(:,32);
  %is = xs(:,[1:31,33:end]);
  is = zeros(length(ts),10);
  for j = 1:length(ts)
      [vain, is(j,:)] = severi_fantini_charawi_difrancesco_2012(ts(j), xs(j,:), p);
  end
  
  peak_time_inds = vs(1:end-2) < vs(2:end-1) & vs(2:end-1) >= vs(3:end) & vs(2:end-1) > PEAK_THR;
  peak_times = ts(peak_time_inds);
  

end
  
function [dy,is] = severi_fantini_charawi_difrancesco_2012(time, states, p)  

  % --- State values --- 

  % --- I cal fca gate ---
  fCa = states(1);

  % --- Ca sr release ---
  I = states(2);
  O = states(3);
  RI = states(4);
  R_Ca = states(5);

  % --- Ca buffering ---
  fCMi = states(6);
  fCMs = states(7);
  fCQ = states(8);
  fTC = states(9);
  fTMC = states(10);
  fTMM = states(11);

  % --- I f y gate ---
  y = states(12);

  % --- I na m gate ---
  m = states(13);

  % --- I na h gate ---
  h = states(14);

  % --- I cal dl gate ---
  dL = states(15);

  % --- I cal fl gate ---
  fL = states(16);

  % --- I cat dt gate ---
  dT = states(17);

  % --- I cat ft gate ---
  fT = states(18);

  % --- Ca dynamics ---
  Ca_jsr = states(19);
  Ca_nsr = states(20);
  Ca_sub = states(21);
  Cai = states(22);
  fBAPTA = states(23);
  fBAPTA_sub = states(24);

  % --- I to q gate ---
  q = states(25);

  % --- I to r gate ---
  r = states(26);

  % --- I kr pa gate ---
  paF = states(27);
  paS = states(28);

  % --- I kr pi gate ---
  piy = states(29);

  % --- I ks n gate ---
  n = states(30);

  % --- I kach a gate ---
  a = states(31);

  % --- Membrane ---
  V_ode = states(32);

  % --- Nai concentration ---
  Nai_ = states(33);

  % Voltage clamp
  V_clamp = (((time > p.t_holding) & (time < p.t_holding +...
    p.t_test))*(p.V_test) + ~((time > p.t_holding) & (time < p.t_holding +...
    p.t_test))*(p.V_holding));

  % I cal fca gate
  fCa_infinity = p.Km_fCa/(p.Km_fCa + Ca_sub);
  tau_fCa = 0.001*fCa_infinity/p.alpha_fCa;

  % Ca sr release
  j_SRCarel = p.ks*O*(Ca_jsr - Ca_sub);
  kCaSR = p.MaxSR - (p.MaxSR - p.MinSR)/(1 + (p.EC50_SR/Ca_jsr)^p.HSR);
  koSRCa = p.koCa/kCaSR;
  kiSRCa = p.kiCa*kCaSR;

  % Ca intracellular fluxes
  b_up = ((p.Iso_1_uM > 0)*(-0.25) + ~(p.Iso_1_uM > 0)*(((p.ACh > 0)*(0.7*p.ACh/(9.0e-5 + p.ACh)) + ~(p.ACh > 0)*(0))));
  P_up = p.P_up_basal*(1 - b_up);
  j_Ca_dif = (Ca_sub - Cai)/p.tau_dif_Ca;
  j_up = P_up/(1 + p.K_up/Cai);
  j_tr = (Ca_nsr - Ca_jsr)/p.tau_tr;

  % Ca buffering
  delta_fTC = p.kf_TC*Cai*(1 - fTC) - p.kb_TC*fTC;
  delta_fTMC = p.kf_TMC*Cai*(1 - fTMC - fTMM) - p.kb_TMC*fTMC;
  delta_fTMM = p.kf_TMM*p.Mgi*(1 - fTMC - fTMM) - p.kb_TMM*fTMM;
  delta_fCMi = p.kf_CM*Cai*(1 - fCMi) - p.kb_CM*fCMi;
  delta_fCMs = p.kf_CM*Ca_sub*(1 - fCMs) - p.kb_CM*fCMs;
  delta_fCQ = p.kf_CQ*Ca_jsr*(1 - fCQ) - p.kb_CQ*fCQ;

  % Cell parameters
  V_cell = 1.0e-9*pi*p.R_cell^2*p.L_cell;
  V_sub = 2.0e-9*pi*p.L_sub*(p.R_cell - p.L_sub/2)*p.L_cell;
  V_jsr = p.V_jsr_part*V_cell;
  V_i = p.V_i_part*V_cell - V_sub;
  V_nsr = p.V_nsr_part*V_cell;

  % Extracted equations
  RTONF = p.R*p.T/p.F;
  V = ((p.clamp_mode >= 1)*(V_clamp) + ~(p.clamp_mode >= 1)*(V_ode));
  Nai = ((p.BAPTA_10_mM > 0)*(7.5) + ~(p.BAPTA_10_mM > 0)*(Nai_));

  % Ionic values
  E_Na = RTONF*log(p.Nao/Nai);
  E_K = RTONF*log(p.Ko/p.Ki);
  E_Ca = 0.5*RTONF*log(p.Cao/Ca_sub);

  % I f
  g_f_Na = ((p.Iva_3_uM >= 1)*(0.03*0.34) + ~(p.Iva_3_uM >= 1)*(0.03));
  g_f_K = ((p.Iva_3_uM >= 1)*(0.03*0.34) + ~(p.Iva_3_uM >= 1)*(0.03));
  ICs_on_Icontrol = ((p.Cs_5_mM >= 1)*(10.6015/5/(10.6015/5 +...
    exp((-0.71*V)/25))) + ~(p.Cs_5_mM >= 1)*(1));
  i_fNa = y^2*p.Ko/(p.Ko + p.Km_f)*g_f_Na*(V - E_Na)*ICs_on_Icontrol;
  i_fK = y^2*p.Ko/(p.Ko + p.Km_f)*g_f_K*(V - E_K)*ICs_on_Icontrol;
  i_f = i_fNa + i_fK;

  % I f y gate
  ACh_shift = ((p.ACh > 0)*(-1 - 9.898*p.ACh^0.618/(p.ACh^0.618 + 0.00122423)) + ~(p.ACh > 0)*(0));
  Iso_shift = ((p.Iso_1_uM > 0)*(7.5) + ~(p.Iso_1_uM > 0)*(0));
  tau_y = p.tauh_f/(0.0708*exp((-V - 5 - (-ACh_shift) - (-Iso_shift))/20.2791) + 10.6*exp(-(ACh_shift + Iso_shift - V)/18));
  y_infinity = 1.0/(1 + exp(-(p.offh_f - V + ACh_shift + Iso_shift)/p.sloh_f));

  % I nak
  Iso_increase = ((p.Iso_1_uM > 0)*(1.2) + ~(p.Iso_1_uM > 0)*(1));
  i_NaK = Iso_increase*p.i_NaK_max/(1 + (p.Km_Kp/p.Ko)^1.2)/(1 +...
    (p.Km_Nap/Nai)^1.3)/(1 + exp((-V - (-E_Na) - 110)/20));

  % I naca
  k43 = Nai/(p.K3ni + Nai);
  k41 = exp((-p.Qn)*V/(2*RTONF));
  di = 1 + Ca_sub/p.Kci*(1 + exp((-p.Qci)*V/RTONF) + Nai/p.Kcni) +...
    Nai/p.K1ni*(1 + Nai/p.K2ni*(1 + Nai/p.K3ni));
  k34 = p.Nao/(p.K3no + p.Nao);
  k32 = exp(p.Qn*V/(2*RTONF));
  do_ = 1 + p.Cao/p.Kco*(1 + exp(p.Qco*V/RTONF)) + p.Nao/p.K1no*(1 +...
    p.Nao/p.K2no*(1 + p.Nao/p.K3no));
  k12 = Ca_sub/p.Kci*exp((-p.Qci)*V/RTONF)/di;
  k14 = Nai/p.K1ni*Nai/p.K2ni*(1 + Nai/p.K3ni)*exp(p.Qn*V/(2*RTONF))/di;
  k21 = p.Cao/p.Kco*exp(p.Qco*V/RTONF)/do_;
  k23 = p.Nao/p.K1no*p.Nao/p.K2no*(1 +...
    p.Nao/p.K3no)*exp((-p.Qn)*V/(2*RTONF))/do_;
  x1 = k41*k34*(k23 + k21) + k21*k32*(k43 + k41);
  x2 = k32*k43*(k14 + k12) + k41*k12*(k34 + k32);
  x3 = k14*k43*(k23 + k21) + k12*k23*(k43 + k41);
  x4 = k23*k34*(k14 + k12) + k14*k21*(k34 + k32);
  i_NaCa = p.K_NaCa*(x2*k21 - x1*k12)/(x1 + x2 + x3 + x4);

  % I na
  E_mh = RTONF*log((p.Nao + 0.12*p.Ko)/(Nai + 0.12*p.Ki));
  i_Na = p.g_Na*m^3*h*(V - E_mh);

  % I na m gate
  if abs(p.offm1_Na-V) < p.delta_m
    alpha_m = 200/p.taum_Na*p.slom1_Na;
  else
    alpha_m = -200/p.taum_Na*(p.offm1_Na-V)/(1 - exp((p.offm1_Na-V)/p.slom1_Na));
  end
  beta_m = 8000/p.taum_Na*exp((p.offm2_Na-V)/p.slom2_Na);

  % I na h gate
  alpha_h = 20/p.tauh_Na*exp((p.offh_Na-V)/p.sloh1_Na);
  beta_h = 2000/p.tauh_Na/(320*exp((p.offh_Na-V)/p.sloh2_Na) + 1);

  % I cal
  Iso_increase = ((p.Iso_1_uM > 0)*(1.23) + ~(p.Iso_1_uM > 0)*(1));
  i_siCa = 2*p.P_CaL*V/(RTONF*(1 - exp((-V)*2/RTONF)))*(Ca_sub -...
    p.Cao*exp((-2*V)/RTONF))*dL*fL*fCa;
  i_siK = 0.000365*p.P_CaL*V/(RTONF*(1 - exp((-V)/RTONF)))*(p.Ki -...
    p.Ko*exp((-V)/RTONF))*dL*fL*fCa;
  i_siNa = 1.85e-5*p.P_CaL*V/(RTONF*(1 - exp((-V)/RTONF)))*(Nai -...
    p.Nao*exp((-V)/RTONF))*dL*fL*fCa;
  ACh_block = 0.31*p.ACh/(p.ACh + 9.0e-5);
  i_CaL = (i_siCa + i_siK + i_siNa)*(1 - ACh_block)*Iso_increase;

  % I cal dl gate
  Iso_shift = ((p.Iso_1_uM > 0)*(-8) + ~(p.Iso_1_uM > 0)*(0));
  Iso_slope = ((p.Iso_1_uM > 0)*(0.69) + ~(p.Iso_1_uM > 0)*(1));
  dL_infinity = 1.0/(1 + exp((p.offm_CaL-V + Iso_shift)/(Iso_slope*p.slom_CaL)));
  adVm = ((V == -41.8)*(-41.80001) + ~(V == -41.8)*(((V == 0)*(0) + ~(V == 0)*(((V == -6.8)*(-6.80001) + ~(V == -6.8)*(V))))));
  bdVm = ((V == -1.8)*(-1.80001) + ~(V == -1.8)*(V));
  alpha_dL = (0.02839*(-41.8-adVm + Iso_shift))/(exp((-41.8-adVm + Iso_shift)/2.5) - 1) + 0.0849*(-6.8-adVm + Iso_shift)/(exp((-6.8-adVm + Iso_shift)/4.8) - 1);
  %alpha_dL = (-0.02839*(adVm + 41.8 - Iso_shift))/(exp((-adVm - 41.8 - (-Iso_shift))/2.5) - 1) - 0.0849*(adVm + 6.8 - Iso_shift)/(exp((-adVm - 6.8 - (-Iso_shift))/4.8) - 1);

  beta_dL = -0.01143*(-1.8-bdVm + Iso_shift)/(exp(-(-1.8-bdVm + Iso_shift)/2.5) - 1);
  %beta_dL = 0.01143*(bdVm + 1.8 - Iso_shift)/(exp((bdVm + 1.8 - Iso_shift)/2.5) - 1);
  tau_dL = p.taum_CaL/(alpha_dL + beta_dL);

  % I cal fl gate
  fL_infinity = 1.0/(1 + exp(-(p.offh_CaL-V)/p.sloh_CaL));
  tau_fL = p.tauh_CaL*(44.3 + 230*exp(-((-36-V)/10)^2));

  % I cat
  i_CaT = 2*p.P_CaT*V/(RTONF*(1 - exp((0-V)*2/RTONF)))*(Ca_sub - p.Cao*exp((-2*V)/RTONF))*dT*fT;

  % I cat dt gate
  dT_infinity = 1.0/(1 + exp((p.offm_CaT-V)/p.slom_CaT));
  tau_dT = p.taum_CaT/(1.068*exp(-(-38.3-V)/30) + 1.068*exp((-38.3-V)/30));

  % I cat ft gate
  fT_infinity = 1.0/(1 + exp(-(p.offh_CaT-V)/p.sloh_CaT));
  tau_fT = p.tauh_CaT/(16.67*exp((-75-V)/83.3) + 16.67*exp(-(-75-V)/15.38));

  % Ca dynamics
  BAPTA = (((p.BAPTA_10_mM > 0) & (time > p.T_Ca))*(10) + ~((p.BAPTA_10_mM > 0) & (time > p.T_Ca))*(0));

  % I to
  i_to = p.g_to*(V - E_K)*q*r;

  % I to q gate
  q_infinity = 1.0/(1 + exp(-(-49-V)/13));
  tau_q = 0.0006*(65.17/(0.57*exp(0.08*(-44-V)) + 0.065*exp(-0.1*(-45.93-V))) + 10.1);

  % I to r gate
  r_infinity = 1.0/(1 + exp((19.3-V)/15));
  tau_r = 0.000924*(15.59/(1.037*exp(-0.09*(-30.61-V)) + 0.369*exp(0.12*(-23.84-V))) + 2.98);

  % I kr
  i_Kr = p.g_Kr*(V - E_K)*(0.9*paF + 0.1*paS)*piy;

  % I kr pa gate
  %(unused!) alfapaF = 1/((1 + exp((-23.2-V)/6.6))*0.84655354)/(37.2*exp(-(0-V)/11.9) + 0.96*exp((0-V)/18.5));
  %(unused!) betapaF = 4*((37.2*exp(-(0-V)/15.9) + 0.96*exp((0-V)/22.5))/0.84655354 - 1/((1 + exp((-23.2-V)/10.6))*0.84655354)/(37.2*exp(-(0-V)/15.9) + 0.96*exp((0-V)/22.5)));
  pa_infinity = 1.0/(1 + exp((p.offm_Kr-V)/p.slom_Kr));
  tau_paS = p.taum1_Kr/(4.2*exp(-(0-V)/17) + 0.15*exp((0-V)/21.6));
  tau_paF = p.taum2_Kr/(30*exp(-(0-V)/10) + exp((0-V)/12));

  % I kr pi gate
  tau_pi = p.tauh_Kr/(100*exp((0-V)/54.645) + 656*exp(-(0-V)/106.157));
  pi_infinity = 1.0/(1 + exp(-(p.offh_Kr-V)/p.sloh_Kr));

  % I ks
  g_Ks = ((p.Iso_1_uM > 0)*(0.00198912) + ~(p.Iso_1_uM > 0)*(0.0016576));
  E_Ks = RTONF*log((p.Ko + 0*p.Nao)/(p.Ki + 0*Nai));
  i_Ks = g_Ks*(V - E_Ks)*n^2;

  % I ks n gate
  Iso_shift = ((p.Iso_1_uM > 0)*(-14) + ~(p.Iso_1_uM > 0)*(0));
  n_infinity = 14/(1 + exp((p.offm1_Ks-V +Iso_shift)/p.slom1_Ks))/(14/(1 + exp((p.offm1_Ks-V +Iso_shift)/p.slom1_Ks)) + exp((Iso_shift+p.offm2_Ks-V)/p.slom2_Ks));
  alpha_n = 28/(1 + exp((40-V+Iso_shift)/3));
  beta_n = exp((5+Iso_shift-V)/25);
  tau_n = p.taum_Ks/(alpha_n + beta_n);

  % I kach
  i_KACh = ((p.ACh > 0)*(p.g_KACh*(V - E_K)*(1 + exp(-(-20-V)/20))*a) +...
    ~(p.ACh > 0)*(0));

  % I kach a gate
  alpha_a = 3.573159/(1 + 1.2155e-6/(p.ACh^1.6951)) + 0.025641;
  beta_a = 10*exp(-0.0133*(-40-V));
  a_infinity = alpha_a/(alpha_a + beta_a);
  tau_a = 1.0/(alpha_a + beta_a);

  % Membrane
  i_tot = i_f + i_Kr + i_Ks + i_to + i_NaK + i_NaCa + i_Na + i_CaL + i_CaT + i_KACh;

  % Nai concentration

  % The ODE system: 33 states

  % Right hand side
  dy = zeros(33, 1);
  dy(1) = 0.001*((fCa_infinity - fCa)/tau_fCa);
  dy(2) = 0.001*(kiSRCa*Ca_sub*O - p.kim*I - p.kom*I - (-koSRCa*Ca_sub^2*RI));
  dy(3) = 0.001*(koSRCa*Ca_sub^2*R_Ca - p.kom*O - kiSRCa*Ca_sub*O -...
    (-p.kim*I));
  dy(4) = 0.001*(p.kom*I - koSRCa*Ca_sub^2*RI - p.kim*RI -...
    (-kiSRCa*Ca_sub*R_Ca));
  dy(5) = 0.001*(p.kim*RI - kiSRCa*Ca_sub*R_Ca - koSRCa*Ca_sub^2*R_Ca -...
    (-p.kom*O));
  dy(6) = 0.001*delta_fCMi;
  dy(7) = 0.001*delta_fCMs;
  dy(8) = 0.001*delta_fCQ;
  dy(9) = 0.001*delta_fTC;
  dy(10) = 0.001*delta_fTMC;
  dy(11) = 0.001*delta_fTMM;
  dy(12) = 0.001*((y_infinity - y)/tau_y);
  dy(13) = 0.001*(alpha_m*(1 - m) - beta_m*m);
  dy(14) = 0.001*(alpha_h*(1 - h) - beta_h*h);
  dy(15) = 0.001*((dL_infinity - dL)/tau_dL);
  dy(16) = 0.001*((fL_infinity - fL)/tau_fL);
  dy(17) = 0.001*((dT_infinity - dT)/tau_dT);
  dy(18) = 0.001*((fT_infinity - fT)/tau_fT);
  dy(19) = 0.001*(j_tr - j_SRCarel - p.CQ_tot*delta_fCQ);
  dy(20) = 0.001*(j_up - j_tr*V_jsr/V_nsr);
  dy(21) = 0.001*(j_SRCarel*V_jsr/V_sub - (i_siCa + i_CaT -...
    2*i_NaCa)/(2*p.F*V_sub) - j_Ca_dif - p.CM_tot*delta_fCMs -...
    p.kfBAPTA*Ca_sub*(BAPTA - fBAPTA_sub) - (-p.kbBAPTA*fBAPTA_sub));
  dy(22) = 0.001*((j_Ca_dif*V_sub - j_up*V_nsr)/V_i - p.CM_tot*delta_fCMi -...
    p.TC_tot*delta_fTC - p.TMC_tot*delta_fTMC - p.kfBAPTA*Cai*(BAPTA -...
    fBAPTA) - (-p.kbBAPTA*fBAPTA));
  dy(23) = 0.001*BAPTA;
  dy(24) = 0.001*BAPTA;
  dy(25) = 0.001*((q_infinity - q)/tau_q);
  dy(26) = 0.001*((r_infinity - r)/tau_r);
  dy(27) = 0.001*((pa_infinity - paF)/tau_paF);
  dy(28) = 0.001*((pa_infinity - paS)/tau_paS);
  dy(29) = 0.001*((pi_infinity - piy)/tau_pi);
  dy(30) = 0.001*((n_infinity - n)/tau_n);
  dy(31) = 0.001*((a_infinity - a)/tau_a);
  dy(32) = 0.001*((-i_tot)/p.C);
  dy(33) = 0.001*((-(i_Na + i_fNa + i_siNa + 3*i_NaK + 3*i_NaCa))/((V_i + V_sub)*p.F));
  is = [i_f, i_Kr, i_Ks, i_to, i_NaK, i_NaCa, i_Na, i_CaL, i_CaT, i_KACh];
end