%
%A mathematical model of action potentials of mouse sinoatrial node cells with molecular bases
%Sanjay Kharche, Jian Yu, Ming Lei, and Henggui Zhang.
%AJP - Heart September 2011 vol. 301 no. 3 H945-H963
%
%LINK TO PAPER:
%
%http://ajpheart.physiology.org/content/301/3/H945.long
%
%Implemented in MATLAB by Tuomo M?ki-Marttunen, 2015

function [vs,is,ts,peak_times] = kharche_SA(T,parChange,isplotted,tol,PEAK_THR,currs)

if nargin < 1 || isempty(T), T = 2000; end
if nargin < 2 || ~isstruct(parChange), if nargin >= 2 && ~isempty(parChange), disp('Warning: parChange not struct, applying no mutations!'); end; parChange = struct; end
if nargin < 3 || isempty(isplotted), isplotted = 0; end
if nargin < 4 || isempty(tol), tol = 1e-9; end
if nargin < 5 || isempty(PEAK_THR), PEAK_THR = 0; end %mV; recognise only local maxima that are above PEAK_THR
if nargin < 6, currs = []; end
if ~isempty(currs) && ~iscell(currs), currs = {currs}; end

R = 8.314472;
Temp = 310.5;
F = 96.4845;

capacitance = 0.025;
vrel = 0.0036;
vsub = 0.03328117;
vup = 0.0348;
vi = 1.34671883;
Mgi = 2.5;
nao = 140.0;
cao = 1.8;
ko = 5.4;
eist = 17.0;
ecal = 47.0;
kmfca = 0.00035;
alpha_fca = 0.021;
ecat = 45.0;
enattxr = 41.5761;
kmnap = 14.0;
kmkp = 1.4;
K1ni = 395.3;
K1no = 1628;
K2ni = 2.289;
K2no = 561.4;
K3ni = 26.44;
K3no = 4.663;
Kci = 0.0207;
Kco = 3.663;
Kcni = 26.44;
Qci = 0.1369;
Qco = 0.0;
Qn = 0.4315;
tdifca = 0.04;

Ttr = 40.0;
ConcTC = 0.031;
ConcTMC = 0.062;
kfTC = 88.8;
kfTMC = 237.7;
kbTC = 0.446;
kbTMC = 0.00751;
kfTMM = 2.277;
kbTMM = 0.751;
ConcCM = 0.045;
kfCM = 237.7;
kbCM = 0.542;
ConcCQ = 10.0;
kfCQ = 0.534;
kbCQ = 0.445;
koca = 10.0;
kom = 0.06;
kica = 0.5;
kim = 0.005;
eca50sr = 0.45;
maxsr = 15.0;
minsr = 1.0;
hsrr = 2.5;
pumphill = 2.0;

FRT = F/(R*Temp);

v = -64.5216286940;
dst = 0.6246780312;
fst = 0.4537033169;
dt = 0.0016256324;
ft = 0.4264459666;
ikr_act = 0.4043600437;
ikr_inact = 0.9250035423;
%ikr_inact2 = 0.1875749806;
iks_act = 0.0127086259;
fl12 = 0.9968141226;
dl12 = 0.0000045583;
fl13 = 0.9809298233;
dl13 = 0.0002036671;
r = 0.0046263658;
m_ttxr = 0.4014088304;
h_ttxr = 0.2724817537;
j_ttxr = 0.0249208708;
m_ttxs = 0.1079085266;
h_ttxs = 0.4500098710;
j_ttxs = 0.0268486392;
y_1_2 = 0.0279984462;
%y_4 = 0.0137659036;
carel = 0.1187281829;
caup = 1.5768287365;
casub = 0.0000560497;
Ftc = 0.0063427103;
Ftmc = 0.1296677919;
Ftmm = 0.7688656371;
Fcms = 0.0242054739;
Fcmi = 0.0138533048;
Fcq = 0.1203184861;
cai = 0.0000319121;
q = 0.6107148187;
fca = 0.7649576191;
nai = 8.1179761505;
ki = 139.8854603066;
resting = 0.7720290515;
open = 0.0000000760;
inactivated = 0.0000000213;
resting_inactivated = 0.2162168926;

Pup_SERCA = 0.04;
ks = 1300000.0;
%Pup2 = 0.04;
pumpkmf = 0.00008;
pumpkmr = 4.5;

ena = (R*Temp/F)*log(nao/nai);
ek  = (R*Temp/F)*log(ko/ki);
eks = ((R*Temp)/F)*log((ko+0.12*nao)/(ki+0.12*nai));
eca = (R*Temp/(2*F))*log(cao/casub);
  

%Channel conductances:
g_st = 0.00006;
g_bNa = 0.0001215;
g_bCa = 0.000015;
g_bK = 0.0000025;
g_K1 = 0.229*0.0039228*0.9;
g_Ks = 0.000299;
g_sus = 0.00039060;
g_NaTTXS = 0.1*5.925e-05;
g_NaTTXR = 0.1*5.925e-05;
g_CaL12 = 0.0010*4.0*1.5;
g_CaL13 = 0.0030*4.0*1.5;
g_CaT = 0.75*0.01862; 
g_f = 0.0057; 
g_Kr = 0.8*0.002955;
g_to = 0.00492;
i_NaK = 0.077;
k_NaCa = 5.5;

%Channel gating parameters:
offm_st = -67.0;
offm_CaT = -26.0;
slom_CaT = 6.0;
taum_CaT = 1.0;
offh_CaT = -61.7;
sloh_CaT = 5.6;
tauh_CaT = 1.0;
offm_Kr = -21.173694;
slom_Kr = 9.757086;
taum_Kr = 0.699821;
offh_Kr = -16.758474;
sloh_Kr = 19.0;
tauh1_Kr = 0.2;
tauh2_Kr = 0.9;
offm_Ks = 20.876040;
slom_Ks = 11.852723;
taum_Ks = 1000.0;
offm_CaL12 = -3.0;
slom_CaL12 = 5.0;
offh_CaL12 = -36;
sloh_CaL12 = 4.6;
offm_CaL13 = -13.5;
slom_CaL13 = 6.0;
offh_CaL13 = -35;
sloh_CaL13 = 7.3;
taum_CaL = 2000;
tauh_CaL = 1.0;
offm_NaTTXS = -31.097331;
slom_NaTTXS = 5.0;
offh_NaTTXS = -56.0;
sloh_NaTTXS = 3.0;
taum_NaTTXS = 1000.0;
tauh_NaTTXS = 1000.0;
tauj_NaTTXS = 1000.0;
offm_NaTTXR = -45.213705;
slom_NaTTXR = 7.219547;
offh_NaTTXR = -62.578120;
sloh_NaTTXR = 6.084036;
taum_NaTTXR = 1000.0;
tauh_NaTTXR = 1000.0;
tauj_NaTTXR = 1000.0;
offh_f = -106.8;
sloh_f = 16.3;
tauh_f = 1.5049;

if isfield(parChange,'offm_st'), offm_st = parChange.offm_st; end
if isfield(parChange,'offm_CaT'), offm_CaT = parChange.offm_CaT; end
if isfield(parChange,'slom_CaT'), slom_CaT = parChange.slom_CaT; end
if isfield(parChange,'taum_CaT'), taum_CaT = parChange.taum_CaT; end
if isfield(parChange,'offh_CaT'), offh_CaT = parChange.offh_CaT; end
if isfield(parChange,'sloh_CaT'), sloh_CaT = parChange.sloh_CaT; end
if isfield(parChange,'tauh_CaT'), tauh_CaT = parChange.tauh_CaT; end
if isfield(parChange,'offm_Kr'), offm_Kr = parChange.offm_Kr; end
if isfield(parChange,'slom_Kr'), slom_Kr = parChange.slom_Kr; end
if isfield(parChange,'taum_Kr'), taum_Kr = parChange.taum_Kr; end
if isfield(parChange,'offh_Kr'), offh_Kr = parChange.offh_Kr; end
if isfield(parChange,'sloh_Kr'), sloh_Kr = parChange.sloh_Kr; end
if isfield(parChange,'tauh1_Kr'), tauh1_Kr = parChange.tauh1_Kr; end
if isfield(parChange,'tauh2_Kr'), tauh2_Kr = parChange.tauh2_Kr; end
if isfield(parChange,'offm_Ks'), offm_Ks = parChange.offm_Ks; end
if isfield(parChange,'slom_Ks'), slom_Ks = parChange.slom_Ks; end
if isfield(parChange,'taum_Ks'), taum_Ks = parChange.taum_Ks; end
if isfield(parChange,'offm_CaL12'), offm_CaL12 = parChange.offm_CaL12; end
if isfield(parChange,'slom_CaL12'), slom_CaL12 = parChange.slom_CaL12; end
if isfield(parChange,'offh_CaL12'), offh_CaL12 = parChange.offh_CaL12; end
if isfield(parChange,'sloh_CaL12'), sloh_CaL12 = parChange.sloh_CaL12; end
if isfield(parChange,'offm_CaL13'), offm_CaL13 = parChange.offm_CaL13; end
if isfield(parChange,'slom_CaL13'), slom_CaL13 = parChange.slom_CaL13; end
if isfield(parChange,'offh_CaL13'), offh_CaL13 = parChange.offh_CaL13; end
if isfield(parChange,'sloh_CaL13'), sloh_CaL13 = parChange.sloh_CaL13; end
if isfield(parChange,'taum_CaL'), taum_CaL = parChange.taum_CaL; end
if isfield(parChange,'tauh_CaL'), tauh_CaL = parChange.tauh_CaL; end
if isfield(parChange,'offm_NaTTXS'), offm_NaTTXS = parChange.offm_NaTTXS; end
if isfield(parChange,'slom_NaTTXS'), slom_NaTTXS = parChange.slom_NaTTXS; end
if isfield(parChange,'offh_NaTTXS'), offh_NaTTXS = parChange.offh_NaTTXS; end
if isfield(parChange,'sloh_NaTTXS'), sloh_NaTTXS = parChange.sloh_NaTTXS; end
if isfield(parChange,'taum_NaTTXS'), taum_NaTTXS = parChange.taum_NaTTXS; end
if isfield(parChange,'tauh_NaTTXS'), tauh_NaTTXS = parChange.tauh_NaTTXS; end
if isfield(parChange,'tauj_NaTTXS'), tauj_NaTTXS = parChange.tauj_NaTTXS; end
if isfield(parChange,'offm_NaTTXR'), offm_NaTTXR = parChange.offm_NaTTXR; end
if isfield(parChange,'slom_NaTTXR'), slom_NaTTXR = parChange.slom_NaTTXR; end
if isfield(parChange,'offh_NaTTXR'), offh_NaTTXR = parChange.offh_NaTTXR; end
if isfield(parChange,'sloh_NaTTXR'), sloh_NaTTXR = parChange.sloh_NaTTXR; end
if isfield(parChange,'taum_NaTTXR'), taum_NaTTXR = parChange.taum_NaTTXR; end
if isfield(parChange,'tauh_NaTTXR'), tauh_NaTTXR = parChange.tauh_NaTTXR; end
if isfield(parChange,'tauj_NaTTXR'), tauj_NaTTXR = parChange.tauj_NaTTXR; end
if isfield(parChange,'offh_f'), offh_f = parChange.offh_f; end
if isfield(parChange,'sloh_f'), sloh_f = parChange.sloh_f; end
if isfield(parChange,'tauh_f'), tauh_f = parChange.tauh_f; end

if isfield(parChange,'g_st'), g_st = parChange.g_st; end
if isfield(parChange,'g_bNa'), g_bNa = parChange.g_bNa; end
if isfield(parChange,'g_bCa'), g_bCa = parChange.g_bCa; end
if isfield(parChange,'g_bK'), g_bK = parChange.g_bK; end
if isfield(parChange,'g_K1'), g_K1 = parChange.g_K1; end
if isfield(parChange,'g_Ks'), g_Ks = parChange.g_Ks; end
if isfield(parChange,'g_sus'), g_sus = parChange.g_sus; end
if isfield(parChange,'g_NaTTXS'), g_NaTTXS = parChange.g_NaTTXS; end
if isfield(parChange,'g_NaTTXR'), g_NaTTXR = parChange.g_NaTTXR; end
if isfield(parChange,'g_CaL12'), g_CaL12 = parChange.g_CaL12; end
if isfield(parChange,'g_CaL13'), g_CaL13 = parChange.g_CaL13; end
if isfield(parChange,'g_CaT'), g_CaT = parChange.g_CaT; end
if isfield(parChange,'g_f'), g_f = parChange.g_f; end
if isfield(parChange,'g_Kr'), g_Kr = parChange.g_Kr; end
if isfield(parChange,'g_to'), g_to = parChange.g_to; end
if isfield(parChange,'i_NaK'), i_NaK = parChange.i_NaK; end
if isfield(parChange,'k_NaCa'), k_NaCa = parChange.k_NaCa; end

if isfield(parChange,'nao'), nao = parChange.nao; end
if isfield(parChange,'cao'), cao = parChange.cao; end
if isfield(parChange,'ko'), ko = parChange.ko; end

inakmax_multiplier = 1.85;
inakmax = inakmax_multiplier*i_NaK;

if isfield(parChange,'Pup_SERCA'), Pup_SERCA = parChange.Pup_SERCA; end
if isfield(parChange,'caup'), caup = parChange.caup; end
if isfield(parChange,'carel'), carel = parChange.carel; end
if isfield(parChange,'vup'), vup = parChange.vup; end
if isfield(parChange,'vrel'), vrel = parChange.vrel; end

data = struct;

data.R = R;
data.Temp = Temp;
data.F = F;
data.capacitance = capacitance;
data.vrel = vrel;
data.vsub = vsub;
data.vup = vup;
data.vi = vi;
data.Mgi = Mgi;
data.nao = nao;
data.cao = cao;
data.ko = ko;
data.eist = eist;
data.ecal = ecal;
data.kmfca = kmfca;
data.alpha_fca = alpha_fca;
data.ecat = ecat;
data.enattxr = enattxr;
data.kmnap = kmnap;
data.kmkp = kmkp;
data.K1ni = K1ni;
data.K1no = K1no;
data.K2ni = K2ni;
data.K2no = K2no;
data.K3ni = K3ni;
data.K3no = K3no;
data.Kci = Kci;
data.Kco = Kco;
data.Kcni = Kcni;
data.Qci = Qci;
data.Qco = Qco;
data.Qn = Qn;
data.tdifca = tdifca;
data.Ttr = Ttr;
data.ConcTC = ConcTC;
data.ConcTMC = ConcTMC;
data.kfTC = kfTC;
data.kfTMC = kfTMC;
data.kbTC = kbTC;
data.kbTMC = kbTMC;
data.kfTMM = kfTMM;
data.kbTMM = kbTMM;
data.ConcCM = ConcCM;
data.kfCM = kfCM;
data.kbCM = kbCM;
data.ConcCQ = ConcCQ;
data.kfCQ = kfCQ;
data.kbCQ = kbCQ;
data.koca = koca;
data.kom = kom;
data.kica = kica;
data.kim = kim;
data.eca50sr = eca50sr;
data.maxsr = maxsr;
data.minsr = minsr;
data.hsrr = hsrr;
data.pumphill = pumphill;
data.FRT = FRT;

data.Pup_SERCA = Pup_SERCA;
data.ks = ks;
data.pumpkmf = pumpkmf;
data.pumpkmr = pumpkmr;
data.g_st = g_st;
data.g_bNa = g_bNa;
data.g_bCa = g_bCa;
data.g_bK = g_bK;
data.g_K1 = g_K1;
data.g_Ks = g_Ks;
data.g_sus = g_sus;
data.g_NaTTXS = g_NaTTXS;
data.g_NaTTXR = g_NaTTXR;
data.g_CaL12 = g_CaL12;
data.g_CaL13 = g_CaL13;
data.g_CaT = g_CaT;
data.g_f = g_f;
data.g_Kr = g_Kr;
data.g_to = g_to;
data.k_NaCa = k_NaCa;
data.offm_st = offm_st;
data.offm_CaT = offm_CaT;
data.slom_CaT = slom_CaT;
data.taum_CaT = taum_CaT;
data.offh_CaT = offh_CaT;
data.sloh_CaT = sloh_CaT;
data.tauh_CaT = tauh_CaT;
data.offm_Kr = offm_Kr;
data.slom_Kr = slom_Kr;
data.taum_Kr = taum_Kr;
data.offh_Kr = offh_Kr;
data.sloh_Kr = sloh_Kr;
data.tauh1_Kr = tauh1_Kr;
data.tauh2_Kr = tauh2_Kr;
data.offm_Ks = offm_Ks;
data.slom_Ks = slom_Ks;
data.taum_Ks = taum_Ks;
data.offm_CaL12 = offm_CaL12;
data.slom_CaL12 = slom_CaL12;
data.offh_CaL12 = offh_CaL12;
data.sloh_CaL12 = sloh_CaL12;
data.offm_CaL13 = offm_CaL13;
data.slom_CaL13 = slom_CaL13;
data.offh_CaL13 = offh_CaL13;
data.sloh_CaL13 = sloh_CaL13;
data.taum_CaL = taum_CaL;
data.tauh_CaL = tauh_CaL;
data.offm_NaTTXS = offm_NaTTXS;
data.slom_NaTTXS = slom_NaTTXS;
data.offh_NaTTXS = offh_NaTTXS;
data.sloh_NaTTXS = sloh_NaTTXS;
data.taum_NaTTXS = taum_NaTTXS;
data.tauh_NaTTXS = tauh_NaTTXS;
data.tauj_NaTTXS = tauj_NaTTXS;
data.offm_NaTTXR = offm_NaTTXR;
data.slom_NaTTXR = slom_NaTTXR;
data.offh_NaTTXR = offh_NaTTXR;
data.sloh_NaTTXR = sloh_NaTTXR;
data.taum_NaTTXR = taum_NaTTXR;
data.tauh_NaTTXR = tauh_NaTTXR;
data.tauj_NaTTXR = tauj_NaTTXR;
data.offh_f = offh_f;
data.sloh_f = sloh_f;
data.tauh_f = tauh_f;
data.inakmax = inakmax;

data.currs = currs;

x = [v,dst,fst,dt,ft,ikr_act,ikr_inact,iks_act,fca,dl13,fl13,dl12,fl12,m_ttxs,h_ttxs,j_ttxs,m_ttxr,h_ttxr,j_ttxr,y_1_2,q,r,...
     resting,open,resting_inactivated,inactivated,Ftc,Ftmc,Ftmm,Fcms,Fcmi,Fcq,casub,cai,carel,caup,nai,ki]';

%dt_int = 0.001; Ndt = T/dt_int; tprev=0; ts = dt_int:dt_int:T; xs = zeros(Ndt,38); is = zeros(Ndt,15); savestep=10;
%for it = 1:Ndt
%  ddt = ts(it) - tprev;
%  [dx, I] = dxdt(ts(it),x,data);
%  x = x + ddt*dx;
%  if mod(it,savestep)==0
%    xs(it/savestep,:) = x';
%    is(it/savestep,:) = I';
%  end  
%  tprev = ts(it);
%end
%ts = ts(savestep:savestep:end);

options = odeset('RelTol',tol,'MaxStep',0.5);
%options = odeset('RelTol',tol);
[ts,xs]=ode15s(@dxdt,[0 T],x,options,data);

vs = xs(:,1);
is = zeros(length(ts),15);
for j = 1:length(ts)
    [vain, is(j,:)] = dxdt(ts(j), xs(j,:), data);
end


if isplotted
  %indss = { {2,3,4,5,7,15}, {1, 6, 11}, {8, 9, 10, 12, 13, 14} };
  %captions = {{'Na-TTX-r','Na-TTX-s','CaL12','CaL13','Kr','to'}, {'f','Ks','CaT'},{'K1', 'st', 'backgr', 'NaK', 'sus', 'Na/Ca'}};
  captions_all = {'f', 'Na-TTX-r','Na-TTX-s','CaL12','CaL13','Ks', 'Kr', 'K1', 'st', 'backgr', 'CaT', 'NaK', 'sus', 'Na/Ca', 'to'};
  indss = { [5 14 11 4 15], [7 12 10 2 13], [1 3 8 9 6] };
  captions = {captions_all(indss{1}), captions_all(indss{2}), captions_all(indss{3})};
  figure;
  plot(t,vs);
  set(gcf,'position',[649   297   202   163])
  figure;
  subplot(4,1,1);
  plot(ts,xs(:,1));ax=axis;axis([t(end)-1000,t(end),ax(3),ax(4)]);
  for iinds=1:3
    inds = indss{iinds};
    subplot(4,1,1+iinds);
    for iind = 1:length(inds)
      disp(num2str(inds(iind)));
      plot(t,sum(xs(:,inds(iind)),2));
      hold on;
    end
    ax=axis;axis([t(end)-1000,t(end),ax(3),ax(4)]);
    legend(captions{iinds});
  end
  set(gcf,'position',[440, 88, 579, 710]);
end

peak_time_inds = xs(1:end-2,1) < xs(2:end-1,1) & xs(2:end-1,1) >= xs(3:end,1) & xs(2:end-1,1) > PEAK_THR;
peak_times = ts(peak_time_inds);



function [dx,I] = dxdt(t,x,data)

  R = data.R;
  Temp = data.Temp;
  F = data.F;
  capacitance = data.capacitance;
  vrel = data.vrel;
  vsub = data.vsub;
  vup = data.vup;
  vi = data.vi;
  Mgi = data.Mgi;
  nao = data.nao;
  cao = data.cao;
  ko = data.ko;
  eist = data.eist;
  ecal = data.ecal;
  kmfca = data.kmfca;
  alpha_fca = data.alpha_fca;
  ecat = data.ecat;
  enattxr = data.enattxr;
  kmnap = data.kmnap;
  kmkp = data.kmkp;
  K1ni = data.K1ni;
  K1no = data.K1no;
  K2ni = data.K2ni;
  K2no = data.K2no;
  K3ni = data.K3ni;
  K3no = data.K3no;
  Kci = data.Kci;
  Kco = data.Kco;
  Kcni = data.Kcni;
  Qci = data.Qci;
  Qco = data.Qco;
  Qn = data.Qn;
  tdifca = data.tdifca;
  Ttr = data.Ttr;
  ConcTC = data.ConcTC;
  ConcTMC = data.ConcTMC;
  kfTC = data.kfTC;
  kfTMC = data.kfTMC;
  kbTC = data.kbTC;
  kbTMC = data.kbTMC;
  kfTMM = data.kfTMM;
  kbTMM = data.kbTMM;
  ConcCM = data.ConcCM;
  kfCM = data.kfCM;
  kbCM = data.kbCM;
  ConcCQ = data.ConcCQ;
  kfCQ = data.kfCQ;
  kbCQ = data.kbCQ;
  koca = data.koca;
  kom = data.kom;
  kica = data.kica;
  kim = data.kim;
  eca50sr = data.eca50sr;
  maxsr = data.maxsr;
  minsr = data.minsr;
  hsrr = data.hsrr;
  pumphill = data.pumphill;
  FRT = data.FRT;

  Pup_SERCA = data.Pup_SERCA;
  ks = data.ks;
  pumpkmf = data.pumpkmf;
  pumpkmr = data.pumpkmr;

  g_st = data.g_st;
  g_bNa = data.g_bNa;
  g_bCa = data.g_bCa;
  g_bK = data.g_bK;
  g_K1 = data.g_K1;
  g_Ks = data.g_Ks;
  g_sus = data.g_sus;
  g_NaTTXS = data.g_NaTTXS;
  g_NaTTXR = data.g_NaTTXR;
  g_CaL12 = data.g_CaL12;
  g_CaL13 = data.g_CaL13;
  g_CaT = data.g_CaT;
  g_f = data.g_f;
  g_Kr = data.g_Kr;
  g_to = data.g_to;
  k_NaCa = data.k_NaCa;
  offm_st = data.offm_st;
  offm_CaT = data.offm_CaT;
  slom_CaT = data.slom_CaT;
  taum_CaT = data.taum_CaT;
  offh_CaT = data.offh_CaT;
  sloh_CaT = data.sloh_CaT;
  tauh_CaT = data.tauh_CaT;
  offm_Kr = data.offm_Kr;
  slom_Kr = data.slom_Kr;
  taum_Kr = data.taum_Kr;
  offh_Kr = data.offh_Kr;
  sloh_Kr = data.sloh_Kr;
  tauh1_Kr = data.tauh1_Kr;
  tauh2_Kr = data.tauh2_Kr;
  offm_Ks = data.offm_Ks;
  slom_Ks = data.slom_Ks;
  taum_Ks = data.taum_Ks;
  offm_CaL12 = data.offm_CaL12;
  slom_CaL12 = data.slom_CaL12;
  offh_CaL12 = data.offh_CaL12;
  sloh_CaL12 = data.sloh_CaL12;
  offm_CaL13 = data.offm_CaL13;
  slom_CaL13 = data.slom_CaL13;
  offh_CaL13 = data.offh_CaL13;
  sloh_CaL13 = data.sloh_CaL13;
  taum_CaL = data.taum_CaL;
  tauh_CaL = data.tauh_CaL;
  offm_NaTTXS = data.offm_NaTTXS;
  slom_NaTTXS = data.slom_NaTTXS;
  offh_NaTTXS = data.offh_NaTTXS;
  sloh_NaTTXS = data.sloh_NaTTXS;
  taum_NaTTXS = data.taum_NaTTXS;
  tauh_NaTTXS = data.tauh_NaTTXS;
  tauj_NaTTXS = data.tauj_NaTTXS;
  offm_NaTTXR = data.offm_NaTTXR;
  slom_NaTTXR = data.slom_NaTTXR;
  offh_NaTTXR = data.offh_NaTTXR;
  sloh_NaTTXR = data.sloh_NaTTXR;
  taum_NaTTXR = data.taum_NaTTXR;
  tauh_NaTTXR = data.tauh_NaTTXR;
  tauj_NaTTXR = data.tauj_NaTTXR;
  offh_f = data.offh_f;
  sloh_f = data.sloh_f;
  tauh_f = data.tauh_f;
  inakmax = data.inakmax;

  currs = data.currs;

  v = x(1);
  dst = x(2);
  fst = x(3);
  dt = x(4);
  ft = x(5);
  ikr_act = x(6);
  ikr_inact = x(7);
  iks_act = x(8);
  fca = x(9);
  dl13 = x(10);
  fl13 = x(11);
  dl12 = x(12);
  fl12 = x(13);
  m_ttxs = x(14);
  h_ttxs = x(15);
  j_ttxs = x(16);
  m_ttxr = x(17);
  h_ttxr = x(18);
  j_ttxr = x(19);
  y_1_2 = x(20);
  q = x(21);
  r = x(22);
  resting = x(23);
  open = x(24);
  resting_inactivated = x(25);
  inactivated = x(26);
  Ftc = x(27);
  Ftmc = x(28);
  Ftmm = x(29);
  Fcms = x(30);
  Fcmi = x(31);
  Fcq = x(32);
  casub = x(33);
  cai = x(34);
  carel = x(35);
  caup = x(36);
  nai = x(37);
  ki = x(38);


  mycurr = 0;
  for istim = 1:length(currs)
    if ~isempty(currs{istim})
      if t >= currs{istim}(1) && t < currs{istim}(2)
        mycurr = mycurr + currs{istim}(3);
      end
    end
  end
  
  %if ceil(rand()*1000)==1
    %disp([num2str(currs{1}) ', t=' num2str(t)])
    %disp([num2str(currs{2}) ', mycurr=' num2str(mycurr)])
  %end
  
  % Ist********************************************************************/

  ena = (R*Temp/F)*log(nao/nai);
  ek  = (R*Temp/F)*log(ko/ki);
  eks = ((R*Temp)/F)*log((ko+0.12*nao)/(ki+0.12*nai));
  eca = (R*Temp/(2*F))*log(cao/casub);
  
  qa = 1.0/(1.0 + exp((offm_st-v)/5.0));
  alphaqa = 1.0/(0.15*exp((0-v)/11.0)+0.2*exp((0-v)/700.0));
  betaqa  =  1.0/(16.0*exp(-(0-v)/8.0)+15.0*exp(-(0-v)/50.0));
  tauqa = 1.0/(alphaqa + betaqa);
  alphaqi = 0.15*1.0/(3100.0*exp(-(-10-v)/13.0)+700.3*exp(-(-10-v)/70.0));
  betaqi =  0.15*1.0/(95.7*exp((-10-v)/10.0) + 50.0*exp((-10-v)/700.0)) + 0.000229/(1+exp((-10-v)/5.0));
  qi = alphaqi/(alphaqi + betaqi);
  tauqi = 1.0/(alphaqi + betaqi);
  ist = g_st*dst*fst*(v - eist); %FIXED: In the original code, ist was calculated basing on the new values of dst and fst, which does not lead to correct numerical integration
  dst_der = ((qa-dst)/tauqa);
  fst_der = ((qi-fst)/tauqi);
  
  % Ib ************************************************************************/
  
  ibna = g_bNa*(v - ena);
  ibca = g_bCa*(v - eca);
  ibk  =  g_bK*(v - ek);
  ib = (ibna + ibca + ibk);
    
  % IK1**********************************************************************/
  
  xk1inf = 1.0/(1.0 + exp(0.070727*(v - ek)));
  ik1 = g_K1*xk1inf*(ko/(ko + 0.228880))*(v - ek);
  
  % ICaT Cav3.1**************************************************************/
  
  tau_dt = taum_CaT/(1.068*exp(-(-26.3-v)/30.0) + 1.068*exp((-26.3-v)/30.0));
  dt_inf = 1.0/(1.0+exp((offm_CaT-v)/slom_CaT));
  tau_ft = tauh_CaT/(0.0153*exp((-61.7-v)/83.3)+0.015*exp(-(-61.7-v)/15.38));
  ft_inf = 1.0/(1.0+exp(-(offh_CaT-v)/sloh_CaT));
  icat = g_CaT*ft*dt*(v - ecat);          
  dt_der = ((dt_inf - dt)/tau_dt);
  ft_der = ((ft_inf - ft)/tau_ft);
  
  % Ikr********************************************************************/
  
  ikr_act_inf = 1.0/(1.0 + exp((offm_Kr-v)/slom_Kr));
  tau_ikr_act = taum_Kr/(0.003596*exp(-(0-v)/15.339290) + 0.000177*exp((0-v)/25.868423));
  
  ikr_inact_inf = 1.0/(1.0 + exp(-(offh_Kr-v)/sloh_Kr));
  tau_ikr_inact = tauh1_Kr+tauh2_Kr/(0.1*exp(-(0-v)/54.645)+0.656*exp(-(0-v)/106.157));
  ikr = g_Kr*ikr_act*ikr_inact*(v - ek);
  ikr_act_der = (ikr_act_inf-ikr_act)/tau_ikr_act;     
  ikr_inact_der = (ikr_inact_inf - ikr_inact)/tau_ikr_inact;
  
  % IKs********************************************************************/
  
  iks_act_inf = 1.0/(1.0 + exp((offm_Ks-v)/slom_Ks));
  tau_iks_act =  taum_Ks/(13.097938/(1.0 + exp((48.910584-v)/10.630272)) + exp((0-v)/35.316539));
  iks = g_Ks*iks_act^2*(v - eks);
  iks_act_der = (iks_act_inf - iks_act)/tau_iks_act;
  
  % ICaL*******************************************************************/
  
  taumcomm_CaL = 1/alpha_fca;
  if abs(v)<=0.001
    alpha_dl  = -28.39*(-(-35.0-v))/(exp((-35.0-v)/2.5)-1.0)+408.173;
  elseif abs(-35.0-v)<=0.001
    alpha_dl  = 70.975+84.9*(0-v)/(exp((0-v)/(1/0.208))-1.0);
  else %if abs(v)>0.001&&fabs(v+35.0)>0.001
    alpha_dl  = 28.39*(-35.0-v)/(exp((-35.0-v)/2.5)-1.0)+84.9*(0-v)/(exp((0-v)/(1/0.208))-1.0);
  end
  
  if abs(5.0-v)<=0.001
   beta_dl   = 28.575;
  else %if fabs(v-5.0)>0.001
    beta_dl   = -11.43*(5.0-v)/(exp(-(5.0-v)/2.5)-1.0);
  end
  
  tau_dl  = taum_CaL/(alpha_dl +beta_dl);
  dl13_inf = 1.0/(1+exp((offm_CaL13-v)/slom_CaL13));
  fl13_inf = 1.0/(1+exp(-(offh_CaL13-v)/sloh_CaL13));
  tau_fl = tauh_CaL*(7.4 + 45.77*exp(-0.5*(-28.1-v)^2/11^2));
  dl12_inf = 1.0/(1+exp((offm_CaL12-v)/slom_CaL12));
  fl12_inf = 1.0/(1+exp(-(offh_CaL12-v)/sloh_CaL12));
  fca_inf = kmfca/(kmfca+casub);
  taufca = fca_inf*taumcomm_CaL;
  ical12 = g_CaL12*fl12*dl12*fca*(v-ecal);
  ical13 = g_CaL13*fl13*dl13*fca*(v-ecal);
  fca_der = (fca_inf - fca)/taufca;
  dl13_der = (dl13_inf - dl13)/tau_dl;
  fl13_der = (fl13_inf - fl13)/tau_fl;
  dl12_der = (dl12_inf - dl12)/tau_dl;
  fl12_der = (fl12_inf - fl12)/tau_fl;
     
  % INa**********************************************************************/
  
  fna = (9.52e-02*exp((-34.4-v)/(1/6.3e-2))/(1+1.66*exp((-63.7-v)/(1/0.225))))+8.69e-2; 
  m3_inf_ttxs = 1.0/(1.0 + exp((offm_NaTTXS-v)/slom_NaTTXS));
  h_inf_ttxs = 1.0/(1.0 + exp(-(offh_NaTTXS-v)/sloh_NaTTXS));
  m_inf_ttxs = m3_inf_ttxs^(1/3); %Using ^(1/3) instead of ^0.333 causes notable effects on beat frequency
  
  tau_m = taum_NaTTXS*((0.6247e-03/(0.832*exp((-56.7-v)/(1/0.335))+0.627*exp(-(-65.01-v)/(1/0.082))))+0.0000492);
  tau_h = tauh_NaTTXS*(((3.717e-06*exp((-17.11-v)/(1/0.2815)))/(1+0.003732*exp((-37.76-v)/(1/0.3426))))+0.0005977);
  tau_j = tauj_NaTTXS*(((0.00000003186*exp((-18.8-v)/(1/0.6219)))/(1+0.00007189*exp((-34.07-v)/(1/0.6683))))+0.003556);
  hs = (1.0-fna)*h_ttxs+fna*j_ttxs;
  
  m3_inf_ttxr = 1.0/(1.0 + exp((offm_NaTTXR-v)/slom_NaTTXR));
  h_inf_ttxr = 1.0/(1.0 + exp(-(offh_NaTTXR-v)/sloh_NaTTXR));
  m_inf_ttxr = m3_inf_ttxr^(1/3); %Using ^(1/3) instead of ^0.333 causes notable effects on beat frequency

  tau_mr = taum_NaTTXR*((0.6247e-03/(0.832*exp((-56.7-v)/(1/0.335))+0.627*exp(-(-65.01-v)/(1/0.082))))+0.0000492);
  tau_hr = tauh_NaTTXR*(((3.717e-06*exp((-17.11-v)/(1/0.2815)))/(1+0.003732*exp((-37.76-v)/(1/0.3426))))+0.0005977);
  tau_jr = tauj_NaTTXR*(((0.00000003186*exp((-18.8-v)/(1/0.6219)))/(1+0.00007189*exp((-34.07-v)/(1/0.6683))))+0.003556);
  hsr = (1.0-fna)*h_ttxr+fna*j_ttxr;
  
  if abs(v)>0.005
    ina_ttxs= g_NaTTXS*m_ttxs^3*hs*nao*(F*F/(R*Temp))*((exp((v-ena)*F/(R*Temp))-1.0)/(exp(v*F/(R*Temp))-1.0))*v;
  else
    ina_ttxs= g_NaTTXS*m_ttxs^3*hs*nao*F*((exp((v-ena)*F/(R*Temp))-1.0));
  end
  
  if abs(v)>0.005
    ina_ttxr = g_NaTTXR*m_ttxr*m_ttxr*m_ttxr*hsr*nao*(F*F/(R*Temp))*((exp((v-enattxr)*F/(R*Temp))-1.0)/(exp(v*F/(R*Temp))-1.0))*v;
  else
    ina_ttxr = g_NaTTXR*m_ttxr*m_ttxr*m_ttxr*hsr*nao*F*((exp((v-enattxr)*F/(R*Temp))-1.0));
  end

  m_ttxs_der = (m_inf_ttxs - m_ttxs)/tau_m;
  h_ttxs_der = (h_inf_ttxs - h_ttxs)/tau_h;
  j_ttxs_der = (h_inf_ttxs - j_ttxs)/tau_j;
  m_ttxr_der = (m_inf_ttxr - m_ttxr)/tau_mr;
  h_ttxr_der = (h_inf_ttxr - h_ttxr)/tau_hr;
  j_ttxr_der = (h_inf_ttxr - j_ttxr)/tau_jr;

  % If**************************************************************************/

  y_inf = 1.0/(1.0 + exp(-(offh_f-v)/sloh_f));
  tau_y_1_2 = tauh_f/(exp((-590.3-v)/(1/0.01094))+ exp(-(85.1-v)/17.2));
  ihk  = 0.6167*g_f*y_1_2*(v - ek);
  ihna = 0.3833*g_f*y_1_2*(v - ena);
  ih = (ihk + ihna);
  y_1_2_der = (y_inf - y_1_2)/tau_y_1_2;
  
  % Ito*************************************************************************/
  q_inf = 1.0/(1.0+exp((v+49.0)/13.0));
  tau_q = (6.06 + 39.102/(0.57*exp((-44.0-v)/12.5)+0.065*exp(-(-45.93-v)/10)))/0.67; 
  r_inf = 1.0/(1.0+exp((19.3-v)/15.0));
  tau_r = (2.75+14.40516/(1.037*exp(-(-30.61-v)/(1/0.09))+0.369*exp((-23.84-v)/(1/0.12))))/0.303;
  ito = g_to*q*r*(v-ek);
  q_der = ((q_inf-q)/tau_q);
  r_der = ((r_inf-r)/tau_r);
   
  % Isus***********************************************************************/
  isus = g_sus*r*(v-ek);
  % Inak***********************************************************************/
  inak = inakmax*(ko^1.2/(kmkp^1.2+ko^1.2))*(nai^1.3/(kmnap^1.3+nai^1.3))/(1.0+exp(-(v-ena+120.0)/30.0));
  
  % iNaCa*******************************************************************/
  
  di=1+(casub/Kci)*(1+exp(-Qci*v*FRT)+nai/Kcni)+(nai/K1ni)*(1+(nai/K2ni)*(1+nai/K3ni));
  doo=1+(cao/Kco)*(1+exp(Qco*v*FRT))+(nao/K1no)*(1+(nao/K2no)*(1+nao/K3no));
  k43=nai/(K3ni+nai);
  k12=(casub/Kci)*exp(-Qci*v*FRT)/di;
  k14=(nai/K1ni)*(nai/K2ni)*(1+nai/K3ni)*exp(Qn*v*FRT/2.0)/di;
  k41=exp(-Qn*v*FRT/2.0);
  k34=nao/(K3no+nao);
  k21=(cao/Kco)*exp(Qco*v*FRT)/doo;
  k23=(nao/K1no)*(nao/K2no)*(1+nao/K3no)*exp(-Qn*v*FRT/2.0)/doo;
  k32=exp(Qn*v*FRT/2);
  x1=k34*k41*(k23+k21)+k21*k32*(k43+k41);
  x2=k43*k32*(k14+k12)+k41*k12*(k34+k32);
  x3=k43*k14*(k23+k21)+k12*k23*(k43+k41);
  x4=k34*k23*(k14+k12)+k21*k14*(k34+k32);
  inaca = k_NaCa*(k21*x2-k12*x1)/(x1+x2+x3+x4);   
  ca_flux = (ical12+ical13+icat-2.0*inaca+ibca)/(2.0*F);
  Jcadif = (casub - cai)/tdifca;
  kcasr = maxsr - (maxsr - minsr)/(1.0 + (eca50sr/carel)^hsrr);
  kosrca = koca/kcasr;
  kisrca = kica*kcasr;
  resting_der             = kim*resting_inactivated - kisrca*casub*resting - kosrca*casub*casub*resting + kom*open;
  open_der                = kosrca*casub*casub*resting - kom*open - kisrca*casub*open + kim*inactivated;
  resting_inactivated_der = kom*inactivated - kosrca*casub*casub*resting_inactivated - kim*resting_inactivated + kisrca*casub*resting;
  inactivated_der         = kisrca*casub*open - kim*inactivated - kom*inactivated + kosrca*casub*casub*resting_inactivated;
  Jrel = ks*open*(carel - casub);
  
  Jup = Pup_SERCA*((cai/pumpkmf)^pumphill - (caup/pumpkmr)^pumphill)/(1.0 + (cai/pumpkmf)^pumphill + (caup/pumpkmr)^pumphill);
  Jtr  = (caup - carel)/Ttr;
  dFtc  = kfTC*cai*(1.0-Ftc)-kbTC*Ftc;
  dFtmc = kfTMC*cai*(1.0-Ftmc-Ftmm)-kbTMC*Ftmc;
  dFtmm = kfTMM*Mgi*(1.0-Ftmc-Ftmm)-kbTMM*Ftmm;
  dFcms = kfCM*casub*(1.0-Fcms)-kbCM*Fcms;
  dFcmi = kfCM*cai*(1.0-Fcmi)-kbCM*Fcmi;
  dFcq  = kfCQ*carel*(1.0-Fcq)-kbCQ*Fcq;
  Ftc_der = dFtc;        
  Ftmc_der = dFtmc;     
  Ftmm_der = dFtmm;     
  Fcms_der = dFcms;     
  Fcmi_der = dFcmi;     
  Fcq_der = dFcq;        
  casub_der = ((-ca_flux+Jrel*vrel)/vsub-Jcadif-ConcCM*dFcms);
  cai_der = ((Jcadif*vsub-Jup*vup)/vi - (ConcCM*dFcmi + ConcTC*dFtc + ConcTMC*dFtmc)); 
  carel_der = (Jtr - Jrel - ConcCQ*dFcq);
  caup_der = (Jup-Jtr*vrel/vup);
  total_current = ih+ina_ttxr+ina_ttxs+ical12+ical13+iks+ikr+ik1+ist+ib+icat+inak+isus+inaca+ito-mycurr;
  v_der = - total_current/capacitance;
  nai_tot = ihna+ina_ttxr+ina_ttxs+3.0*inak+3.0*inaca+ist+ibna;
  ki_tot = ihk+iks+ikr+ik1+ibk-2.0*inak+isus+ito;
  nai_der = -(nai_tot)/(F*vi);
  ki_der = -(ki_tot)/(F*vi);

  dx = zeros(38,1);
  dx(1) = v_der;
  dx(2) = dst_der;
  dx(3) = fst_der;
  dx(4) = dt_der;
  dx(5) = ft_der;
  dx(6) = ikr_act_der;
  dx(7) = ikr_inact_der;
  dx(8) = iks_act_der;
  dx(9) = fca_der;
  dx(10) = dl13_der;
  dx(11) = fl13_der;
  dx(12) = dl12_der;
  dx(13) = fl12_der;
  dx(14) = m_ttxs_der;
  dx(15) = h_ttxs_der;
  dx(16) = j_ttxs_der;
  dx(17) = m_ttxr_der;
  dx(18) = h_ttxr_der;
  dx(19) = j_ttxr_der;
  dx(20) = y_1_2_der;
  dx(21) = q_der;
  dx(22) = r_der;
  dx(23) = resting_der;
  dx(24) = open_der;
  dx(25) = resting_inactivated_der;
  dx(26) = inactivated_der;
  dx(27) = Ftc_der;
  dx(28) = Ftmc_der;
  dx(29) = Ftmm_der;
  dx(30) = Fcms_der;
  dx(31) = Fcmi_der;
  dx(32) = Fcq_der;
  dx(33) = casub_der;
  dx(34) = cai_der;
  dx(35) = carel_der;
  dx(36) = caup_der;
  dx(37) = nai_der;
  dx(38) = ki_der;

  I = [ih, ina_ttxr, ina_ttxs, ical12, ical13, iks, ikr, ik1, ist, ib, icat, inak, isus, inaca, ito];
