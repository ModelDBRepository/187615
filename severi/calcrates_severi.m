function [nspikesAll, nspikes_control] = calcrates_severi(considered,scalingfile)

T = 10000;

tobesaved = 0;

if nargin < 1 || isempty(considered), considered = [1 48 66 78 83 87 91 93 96]; tobesaved=1; end
if nargin < 2 || isempty(scalingfile), scalingfile = 'scalings_severi.mat'; end


MT = getMT_severi();
defVals = getDefVals_severi();
geneNames = getgenenames();

control_file = load('severi_control.mat');
vs_control = control_file.vs';
ts_control = control_file.ts;
lastBelow3000 = max(find(ts_control<3000));
vs3000_control = vs_control(1:lastBelow3000+1);
ts3000_control = ts_control(1:lastBelow3000+1);

dvs_control=membpotderivs(ts_control,vs_control);
Nextrasteps = 21;
peak_times_last = control_file.peak_times(end-1:end);
if peak_times_last(2) > T-0.1
  peak_times_last = control_file.peak_times(end-2:end-1);
end
peak_times_last_inds = zeros(1,2);
for i=1:2, [~,peak_times_last_inds(i)] = min(abs(ts_control-peak_times_last(i))); end
vs_control_lc = vs_control(peak_times_last_inds(1)+1:peak_times_last_inds(2)-1+Nextrasteps);
dvs_control_lc = dvs_control(peak_times_last_inds(1):peak_times_last_inds(2)-2+Nextrasteps);
%vs_control_lc = vs_control_lc(10:10:end); dvs_control_lc = dvs_control_lc(10:10:end);

cols = [0.2667, 0.2667, 0.2667; 0.0039, 0.1373, 0.2706; 0.6667, 0, 0.6667; 0.7333, 0.6667, 0.6667; 0.9333, 0.4, 0;1, 0, 0; 0, 0.6, 0.6;0.4667, 0.1333, 0.4667; 0, 0.8, 0];
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
    nspikesMut = [];
    for iallmutval = 1:cumprodnVals(end)
      for imutvar = 1:length(thisMut)
        if imutvar==1
          allmutvals{iallmutval}(imutvar) = thisMut{imutvar}{2}{1}(mod(iallmutval-1,nVals(imutvar))+1);
        else
          allmutvals{iallmutval}(imutvar) = thisMut{imutvar}{2}{1}(mod(floor((iallmutval-1)/cumprodnVals(imutvar-1)),nVals(imutvar))+1);          
        end
      end
    end
    for iallmutval = 1:cumprodnVals(end)
      counter = counter + 1;
      
      iters = [1,3,6,7,9];
      nspikes = nan(size(iters));
      parChangeRel = struct;
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
              if allmutvals{iallmutval}(imutvar) > 0 && iter==6
                parChangeRel = setfield(parChangeRel,thesemutvars{imutvar}{kmutvar},['+' num2str(thisCoeff*allmutvals{iallmutval}(imutvar))]);
              elseif iter==6
                parChangeRel = setfield(parChangeRel,thesemutvars{imutvar}{kmutvar},[num2str(thisCoeff*allmutvals{iallmutval}(imutvar))]);
              end
            end
          else
            for kmutvar = 1:length(thesemutvars{imutvar})
              parChange = setfield(parChange,thesemutvars{imutvar}{kmutvar},getfield(defVals,thesemutvars{imutvar}{kmutvar}) * allmutvals{iallmutval}(imutvar)^thisCoeff);
              if iter==6
                parChangeRel = setfield(parChangeRel,thesemutvars{imutvar}{kmutvar},['\times' num2str(allmutvals{iallmutval}(imutvar)^thisCoeff)]);
              end
            end        
          end
        end
        parChange
        if ~exist(['severi_' num2str(igene) '_' num2str(imut) '_' num2str(iallmutval) '_' num2str(iter) '.mat'])
          %continue;
          [vs,is,ts,peak_times] = severi_SA(T,parChange);
          vs0 = vs';
          ts0 = ts;
          ts = 0.1:0.1:T;
          vs = interpolate_multidim(ts0,vs0,ts);
          dvs=membpotderivs(ts,vs);          

          lastBelow3000 = max(find(ts<3000));
          vs3000 = vs(1:lastBelow3000+1); ts3000 = ts(1:lastBelow3000+1);

          Nextrasteps = 21;
          
          vs_lc = [];
          dvs_lc = [];
          if length(peak_times) > 1
            peak_times_last = peak_times(end-1:end);
            if peak_times_last(2) > T-0.1
              peak_times_last = peak_times(end-2:end-1);
            end
            peak_times_last_inds = zeros(1,2);
            for i=1:2, [~,peak_times_last_inds(i)] = min(abs(ts-peak_times_last(i))); end
            vs_lc = vs(peak_times_last_inds(1)+1:min(peak_times_last_inds(2)-1+Nextrasteps,length(vs)));
            dvs_lc = dvs(peak_times_last_inds(1):min(peak_times_last_inds(2)-2+Nextrasteps,length(dvs)));
            %vs_lc = vs_lc(10:10:end); dvs_lc = dvs_lc(10:10:end);
            if length(dvs_lc) < length(vs_lc)
              vs_lc = vs_lc(1:length(dvs_lc));
            end
          end
            
          if isempty(vs_lc) || isempty(dvs_lc)
            disp('Empty vs_lc!');
            vs_lc = mean(vs_control_lc); dvs_lc = 0;
          end
          save(['severi_' num2str(igene) '_' num2str(imut) '_' num2str(iallmutval) '_' num2str(iter) '.mat'],'vs3000','ts3000','vs_lc','dvs_lc','peak_times');
        else
          load(['severi_' num2str(igene) '_' num2str(imut) '_' num2str(iallmutval) '_' num2str(iter) '.mat'],'vs3000','ts3000','vs_lc','dvs_lc','peak_times');
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

if ~exist(['severi_control.mat'])
    disp(['severi_control.mat not found'])
else
    load(['severi_control.mat'],'peak_times');
end
        
disp(num2str(peak_times));
npeaks = length(peak_times);
peakfraction = (T-peak_times(end))/(peak_times(end)-peak_times(end-1));
npeaks = npeaks + peakfraction;
nspikes_control = npeaks;

if tobesaved
  save('rates_severi.mat','nspikesAll','nspikes_control');
end
