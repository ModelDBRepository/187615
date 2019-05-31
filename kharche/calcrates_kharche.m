function [nspikesAll,nspikes_control]=calcrates_kharche(considered,scalingfile)

T = 10000;
dt_int = 0.0025;
dt_save = 0.05;
savestep = 4;
skip = round(dt_save/dt_int/savestep);
tobesaved = 0;

if nargin < 1 || isempty(considered), considered = [1 48 66 78 83 87 91 93 96]; tobesaved=1; end
if nargin < 2 || isempty(scalingfile), scalingfile = 'scalings_kharche.mat'; end

MT = getMT_kharche();
defVals = getDefVals_kharche();

control = load('kharche_control.mat');
vs_control = control.vs';
ts_control = control.ts;
dvs_control=membpotderivs(ts_control,vs_control');
Nextrasteps = 21;
peak_times_last = control.peak_times(end-1:end);
if peak_times_last(2) > T-0.1
  peak_times_last = control.peak_times(end-2:end-1);
end
peak_times_last_inds = zeros(1,2);
for i=1:2, [~,peak_times_last_inds(i)] = min(abs(ts_control-peak_times_last(i))); end
vs_control_lc = vs_control(peak_times_last_inds(1)+1:peak_times_last_inds(2)-1+Nextrasteps);
dvs_control_lc = dvs_control(peak_times_last_inds(1):peak_times_last_inds(2)-2+Nextrasteps);
vs_control_lc = vs_control_lc(10:10:end); dvs_control_lc = dvs_control_lc(10:10:end);

cols = [0.2667, 0.2667, 0.2667; 0.0039, 0.1373, 0.2706; 0.6667, 0.6667, 0.6667; 0.7333, 0.6667, 0.6667; 0.9333, 0.4, 0;1, 0, 0; 0, 0.6, 0.6;0.4667, 0.1333, 0.4667; 0, 0.8, 0];
%['#444444','#012345','#aa00aa','#bbaa00','#ee6600','#ff0000', '#009999','#772277','#00cc00']
col_control = [0.1333, 0.1333, 1];
coeffCoeffs = {[0.25,0],[0.125,0],[0.5,0],[0.5,1.0/3],[0.5,2.0/3],[0.5,1.0],[-0.25,0],[-0.125,0],[-0.5,0]};

scaling = load(scalingfile);
theseCoeffsAll = scaling.coeffs(1,:);

counter = 0;
nspikesAll = [];
for igene = 1:length(MT)
  nspikesGene = [];
  for imut = 1:length(MT{igene})
    thisMut = MT{igene}{imut};
    nVals = cellfun(@(x)length(x{2}),thisMut);
    cumprodnVals = cumprod(nVals);
    disp(cumprodnVals)
    thesemutvars = [];
    
    if iscell(theseCoeffsAll{igene})
      theseCoeffs = theseCoeffsAll{igene}{imut};
    else
      if size(theseCoeffsAll{igene},1) > size(theseCoeffsAll{igene},2)
        theseCoeffs = theseCoeffsAll{igene}(imut);
      else
        theseCoeffs = theseCoeffsAll{igene};
      end
    end
    for imutvar=1:length(thisMut)
      if ~iscell(thisMut{imutvar}{1})
        thisMut{imutvar}{1} = {thisMut{imutvar}{1}};
      end
      if ~iscell(thisMut{imutvar}{2})
        thisMut{imutvar}{2} = {thisMut{imutvar}{2}};
      end
      thesemutvars = [thesemutvars, {thisMut{imutvar}{1}}];
    end

    allmutvars = repmat({thesemutvars},cumprodnVals(end),1)
    allmutvals = [];
    for iallmutval = 1:cumprodnVals(end)
      allmutvals = [allmutvals, {zeros(length(thesemutvars),1)}];
    end
    for iallmutval = 1:cumprodnVals(end)
      for imutvar = 1:length(thisMut)
        if imutvar==1
          allmutvals{iallmutval}(imutvar) = thisMut{imutvar}{2}{1}(mod(iallmutval-1,nVals(imutvar))+1);
        else
          allmutvals{iallmutval}(imutvar) = thisMut{imutvar}{2}{1}(mod(floor((iallmutval-1)/cumprodnVals(imutvar-1)),nVals(imutvar))+1);          
        end
      end
    end
    nspikesMut = [];
    for iallmutval = 1:cumprodnVals(end)
      counter = counter + 1;

      iters = [1,3,6,7,9];
      nspikes = nan(size(iters));
      for iiter = 1:length(iters)
        if all(counter ~= considered)
          continue
        end
      
        iter = iters(iiter);
        thisCoeff = coeffCoeffs{iter}(1)*theseCoeffs(iallmutval) + coeffCoeffs{iter}(2)*(1.0 - 0.5*theseCoeffs(iallmutval))
        
        parChange = struct;
        for imutvar = 1:length(thesemutvars)
          if length(strfind(thesemutvars{imutvar}{1},'off')) > 0
            for kmutvar = 1:length(thesemutvars{imutvar})
              parChange = setfield(parChange,thesemutvars{imutvar}{kmutvar},getfield(defVals,thesemutvars{imutvar}{kmutvar}) + thisCoeff*allmutvals{iallmutval}(imutvar));
            end
          else
            for kmutvar = 1:length(thesemutvars{imutvar})
              parChange = setfield(parChange,thesemutvars{imutvar}{kmutvar},getfield(defVals,thesemutvars{imutvar}{kmutvar}) * allmutvals{iallmutval}(imutvar)^thisCoeff);
            end        
          end
        end
        parChange
        if ~exist(['kharche_' num2str(igene) '_' num2str(imut) '_' num2str(iallmutval) '_' num2str(iter) '.mat'])
          disp(['kharche_' num2str(igene) '_' num2str(imut) '_' num2str(iallmutval) '_' num2str(iter) '.mat not found'])
          [vs,is,ts,peak_times] = kharche_SA(T,parChange,0,1e-12);
          dvs=membpotderivs(ts,vs);     
          ts3000 = 0.01:0.01:3000;
          vs3000 = interpolate(ts,vs,ts3000);

          Nextrasteps = 21;
          vs_lc = [];
          dvs_lc = [];
          if length(peak_times) > 1
            peak_time_inds = [length(peak_times)-1, length(peak_times)];
            peak_times_last = peak_times(peak_time_inds);                                                                                                                                                    
            while peak_times_last(2) > ts(end-Nextrasteps-1)                                                                                                                                                                 
              peak_time_inds = peak_time_inds - 1;                                                                                                                                                           
              peak_times_last = peak_times(peak_time_inds);                                                                                                                                                  
            end
            peak_times_last_inds = zeros(1,2);
            for i=1:2, [~,peak_times_last_inds(i)] = min(abs(ts-peak_times_last(i))); end
            vs_lc = vs(peak_times_last_inds(1)+1:peak_times_last_inds(2)-1+Nextrasteps);
            dvs_lc = dvs(peak_times_last_inds(1):peak_times_last_inds(2)-2+Nextrasteps);
            vs_lc = vs_lc(10:10:end); dvs_lc = dvs_lc(10:10:end);
          end
          if isempty(vs_lc) || isempty(dvs_lc)
            disp('Empty vs_lc!');
            vs_lc = mean(vs_control_lc); dvs_lc = 0;
          end
          save(['kharche_' num2str(igene) '_' num2str(imut) '_' num2str(iallmutval) '_' num2str(iter) '.mat'],'vs3000','ts3000','vs_lc','dvs_lc','peak_times');
        else
          load(['kharche_' num2str(igene) '_' num2str(imut) '_' num2str(iallmutval) '_' num2str(iter) '.mat'],'peak_times');
        end
        
        disp(num2str(peak_times));
        npeaks = length(peak_times);
        if npeaks > 6
          peakfraction = (T-peak_times(end))/(peak_times(end)-peak_times(end-1));
          npeaks = npeaks + peakfraction;
        end
        nspikes(iiter) = npeaks;
      end
      nspikesMut = [nspikesMut, {nspikes}];
    end
    nspikesGene = [nspikesGene, {nspikesMut}];
  end
  nspikesAll = [nspikesAll, {nspikesGene}];
end

if ~exist(['kharche_control.mat'])
    disp(['kharche_control.mat not found'])
else
    load(['kharche_control.mat'],'peak_times');
end
        
disp(num2str(peak_times));
npeaks = length(peak_times);
peakfraction = (T-peak_times(end))/(peak_times(end)-peak_times(end-1));
npeaks = npeaks + peakfraction;
nspikes_control = npeaks;

if tobesaved
  save('rates_kharche.mat','nspikesAll','nspikes_control');
end
