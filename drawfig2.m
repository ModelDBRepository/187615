variants = {[0,0,0],[2,3,0],[1,2,13],[3,0,1],[5,0,0]};
variants_matlab = {[1,1,1],[3,4,1],[2,3,14],[4,1,2],[6,1,1]};
variant_chosen = 1;

geneNames = getgenenames();
titles = {'CACNA1C [27]', 'CACNA1D [36]', 'CACNB2 [41]', 'CACNA1I [43]', 'ATP2A2 [45, 46]'};
MT = getMT_kharche();

SP = subplottight2(4,1+length(variants));
SP = SP(end:-1:1,:);

SPinset = zeros(2,length(variants));
for iy=1:4
  for ix=1:length(variants)
    if iy < 3
      SPinset(iy,ix) = axes('position',[0.39+(ix-1)*0.135,0.79-(iy-1)*0.2,0.045,0.07]);
    end
    set(SP(iy,1+ix),'position', [0.31+(ix-1)*0.135, 0.7725-(iy-1)*0.2, 0.135, 0.16])
  end
  set(SP(iy,1),'position', [0.09, 0.7725-(iy-1)*0.2, 0.16, 0.16])
end

data_kharche = struct; data_severi = struct; data_hay = struct; data_almog = struct;
if exist('kharche/rates_kharche.mat')
  data_kharche = load('kharche/rates_kharche.mat'); data_kharche.times=[]; data_kharche.Vsoma=[];
end
if exist('severi/rates_severi.mat')
  data_severi = load('severi/rates_severi.mat'); data_severi.times=[]; data_severi.Vsoma=[];
end
if exist('hay/fig1_curves.mat')
  data_hay = load('hay/fig2_curves.mat');
end
if exist('almog/fig1_curves.mat')
  data_almog = load('almog/fig2_curves.mat');
end

i=variant_chosen;
for iter=[1,3,7,9]
  if exist(['kharche/kharche_' num2str(variants_matlab{i}(1)) '_' num2str(variants_matlab{i}(2)) '_' num2str(variants_matlab{i}(3)) '_' num2str(iter) '.mat'])
    A=load(['kharche/kharche_' num2str(variants_matlab{i}(1)) '_' num2str(variants_matlab{i}(2)) '_' num2str(variants_matlab{i}(3)) '_' num2str(iter) '.mat']);
    tinds = find(A.ts3000 >= 0 & A.ts3000 <= 600);
    data_kharche.times = [data_kharche.times; A.ts3000(tinds)];
    data_kharche.Vsoma = [data_kharche.Vsoma; A.vs3000(tinds)];
  else
    disp(['kharche/kharche_' num2str(variants_matlab{i}(1)) '_' num2str(variants_matlab{i}(2)) '_' num2str(variants_matlab{i}(3)) '_' num2str(iter) '.mat does not exist'])
  end
  A=load(['kharche/kharche_control.mat']);
  tinds = find(A.ts >= 0 & A.ts <= 600);
  data_kharche.times_control = A.ts(tinds);
  data_kharche.Vsoma_control = A.vs(tinds);

  if exist(['severi/severi_' num2str(variants_matlab{i}(1)) '_' num2str(variants_matlab{i}(2)) '_' num2str(variants_matlab{i}(3)) '_' num2str(iter) '.mat'])
    A=load(['severi/severi_' num2str(variants_matlab{i}(1)) '_' num2str(variants_matlab{i}(2)) '_' num2str(variants_matlab{i}(3)) '_' num2str(iter) '.mat']);
    tinds = find(A.ts3000 >= 0 & A.ts3000 <= 600);
    data_severi.times = [data_severi.times; A.ts3000(tinds)];
    data_severi.Vsoma = [data_severi.Vsoma; A.vs3000(tinds)];
  else
    disp( ['severi/severi_' num2str(variants_matlab{i}(1)) '_' num2str(variants_matlab{i}(2)) '_' num2str(variants_matlab{i}(3)) '_' num2str(iter) '.mat does not exist'])
  end
  A=load(['severi/severi_control.mat']);
  tinds = find(A.ts >= 0 & A.ts <= 600);
  data_severi.times_control = A.ts(tinds);
  data_severi.Vsoma_control = A.vs(tinds);
end


cols = [0.2667, 0.2667, 0.2667; 0.0039, 0.1373, 0.2706; 0.6667, 0, 0.6667; 0.7333, 0.6667, 0; 0.9333, 0.4, 0;1, 0, 0; 0, 0.6, 0.6;0.4667, 0.1333, 0.4667; 0, 0.8, 0];
cols = cols([1,3,7,9],:);
col_control = [0.1333, 0.1333, 1];

xs_bars = [2,1,4,5]; x_bar_control = 3;
for iiter=1:4
  iiterplus = iiter+(iiter>2);
  try
    subplot(SP(1,1)); plot(data_hay.times{variant_chosen,iiter},data_hay.Vsoma{variant_chosen,iiter},'color',cols(iiter,:)); hold on;
  catch
    disp('Hay model time course data missing')
  end
  try
    subplot(SP(2,1)); plot(data_almog.times{variant_chosen,iiter},data_almog.Vsoma{variant_chosen,iiter},'color',cols(iiter,:)); hold on;
  catch
    disp('Almog model time course data missing')
  end
  try
    subplot(SP(3,1)); plot(data_kharche.times(iiter,:),data_kharche.Vsoma(iiter,:),'color',cols(iiter,:)); hold on;
  catch
    disp('Kharche model time course data missing')
  end
  try
    subplot(SP(4,1)); plot(data_severi.times(iiter,:),data_severi.Vsoma(iiter,:),'color',cols(iiter,:)); hold on;
  catch
    disp('Severi model time course data missing')
  end
  
  for ivar=1:length(variants)
    try
      subplot(SP(1,1+ivar)); plot(data_hay.Is,reshape(data_hay.spikeFreqs(ivar,iiter,:),size(data_hay.Is)),'color',cols(iiter,:)); hold on;
    catch
      disp('Hay model f-I data missing')
    end    
    try
      subplot(SP(2,1+ivar)); plot(data_almog.Is,reshape(data_almog.spikeFreqs(ivar,iiter,:),size(data_almog.Is)),'color',cols(iiter,:)); hold on;
    catch
      disp('Almog model f-I data missing')
    end 
    try
      subplot(SPinset(1,ivar)); a = bar(xs_bars(iiter),(data_hay.threshIs(ivar,iiter)-data_hay.threshI_control)/data_hay.threshI_control*100); set(a,'facecolor',cols(iiter,:)); hold on;
    catch
      disp('Hay model threshold current data missing')
    end
    try
      subplot(SPinset(2,ivar)); a = bar(xs_bars(iiter),(data_almog.threshIs(ivar,iiter)-data_almog.threshI_control)/data_almog.threshI_control*100); set(a,'facecolor',cols(iiter,:)); hold on;
    catch
      disp('Almog model threshold current data missing')
    end 
    try
      subplot(SP(3,1+ivar)); a = bar(xs_bars(iiter), (data_kharche.nspikesAll{variants_matlab{ivar}(1)}{variants_matlab{ivar}(2)}{variants_matlab{ivar}(3)}(iiterplus)-data_kharche.nspikes_control)/data_kharche.nspikes_control*100); set(a,'facecolor',cols(iiter,:)); hold on;%axis([0,600,-90,40]);
    catch
      disp('Kharche rates data missing')
    end 
    try
      subplot(SP(4,1+ivar)); a = bar(xs_bars(iiter), (data_severi.nspikesAll{variants_matlab{ivar}(1)}{variants_matlab{ivar}(2)}{variants_matlab{ivar}(3)}(iiterplus)-data_severi.nspikes_control)/data_severi.nspikes_control*100); set(a,'facecolor',cols(iiter,:)); hold on;%axis([0,600,-90,40]);
    catch
      disp('Severi rates data missing')
    end 
  end
end

subplot(SP(1,1)); plot(data_hay.times_control,data_hay.Vsoma_control,'color',col_control); hold on; axis([0,600,-90,40]); 
title(titles{variant_chosen});
set(gca,'ytick',-80:40:40); ylabel('mV');
text(-260,-43,'Hay','rot',90,'fontweight','bold');
subplot(SP(2,1)); plot(data_almog.times_control,data_almog.Vsoma_control,'color',col_control); hold on; axis([0,600,-90,40]); set(gca,'ytick',-80:40:40); ylabel('mV');
text(-260,-55,'Almog','rot',90,'fontweight','bold');
subplot(SP(3,1)); plot(data_kharche.times_control,data_kharche.Vsoma_control,'color',col_control); hold on;axis([0,600,-90,40]); set(gca,'ytick',-80:40:40); ylabel('mV');
text(-260,-60,'Kharche','rot',90,'fontweight','bold');
subplot(SP(4,1)); plot(data_severi.times_control,data_severi.Vsoma_control,'color',col_control); hold on;axis([0,600,-90,40]); set(gca,'ytick',-80:40:40); xlabel('ms'); ylabel('mV');
text(-260,-57,'Severi','rot',90,'fontweight','bold');

for ivar=1:length(variants)
  try
    subplot(SP(1,1+ivar)); plot(data_hay.Is,reshape(data_hay.spikeFreqs_control,size(data_hay.Is)),'color',col_control); axis([0.25,1.4,0,20]); set(gca,'xtick',[0.4,0.8,1.2]);
  catch
    disp('Hay model f-I data missing')
  end    
  title(titles{ivar});
  if ivar==1, ylabel('Hz   ','rot',0); end
  try
    subplot(SP(2,1+ivar)); plot(data_almog.Is,reshape(data_almog.spikeFreqs_control,size(data_almog.Is)),'color',col_control); axis([0.7,0.9,0,30]); set(gca,'xtick',[0.7,0.75,0.8,0.85]); if ivar==1, ylabel('Hz   ','rot',0); end
  catch
    disp('Almog model f-I data missing')
  end    
  subplot(SPinset(1,ivar)); a = bar(x_bar_control,0); set(a,'facecolor',col_control); axis([0,6,-1.5,1.7]); set(gca,'xtick',[],'fontsize',5); ylabel('%','rot',0);
  subplot(SPinset(2,ivar)); a = bar(x_bar_control,0); set(a,'facecolor',col_control); axis([0,6,-1.5,1.5]); set(gca,'xtick',[],'fontsize',5); ylabel('%','rot',0);
  subplot(SP(3,1+ivar)); a = bar(x_bar_control, 0); set(a,'facecolor',col_control); axis([0,6,-10,10]); set(gca,'xtick',[]); if ivar==1, ylabel('%','rot',0); end
  subplot(SP(4,1+ivar)); a = bar(x_bar_control, 0); set(a,'facecolor',col_control); axis([0,6,-10,10]); set(gca,'xtick',[]); if ivar==1, ylabel('%','rot',0); end

  mystr = MTentrystr(MT{variants_matlab{ivar}(1)}{variants_matlab{ivar}(2)},variants_matlab{ivar}(3),{'*{\it{c}}','^{\it{c}}'});
  mystr = strrep(mystr,'*','\times');
  t = annotation(gcf,'textbox',[0.31+(ivar-1)*0.135,0.02,0.18,0.14]);set(t,'string',mystr,'fontsize',9,'linestyle','none');
  if ivar>1
    for imod=1:4
      set(SP(imod,1+ivar), 'yticklabel', []);
    end
  end
end    
mystr = MTentrystr(MT{variants_matlab{variant_chosen}(1)}{variants_matlab{variant_chosen}(2)},variants_matlab{variant_chosen}(3),{'*{\it{c}}','^{\it{c}}'});
mystr = strrep(mystr,'*','\times');

t = annotation(gcf,'textbox',[0.08,0.0,0.18,0.12]);set(t,'string',mystr,'fontsize',9,'linestyle','none');
t = annotation(gcf,'textbox',[0.01,0.97,0.1,0.02]);set(t,'string','A','fontsize',24,'linestyle','none');
t = annotation(gcf,'textbox',[0.24,0.97,0.1,0.02]);set(t,'string','B','fontsize',24,'linestyle','none');



print('-depsc','fig2')
