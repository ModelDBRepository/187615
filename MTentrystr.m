function str=MTentrystr(entry,myallmutval,coeff)

if nargin < 2
    myallmutval = 1;
end
if nargin < 3
    coeff = {'',''};
end

nVals = cellfun(@(x)length(x{2}),entry);
cumprodnVals = cumprod(nVals);
thesemutvars = [];
for imutvar=1:length(entry)
  if ~iscell(entry{imutvar}{1})
    entry{imutvar}{1} = {entry{imutvar}{1}};
  end
  if ~iscell(entry{imutvar}{2})
    entry{imutvar}{2} = {entry{imutvar}{2}};
  end
  thesemutvars = [thesemutvars {entry{imutvar}{1}}];
end

allmutvars = repmat({thesemutvars},cumprodnVals(end),1);
allmutvals = [];
for iallmutval = 1:cumprodnVals(end)
  allmutvals = [allmutvals, {zeros(length(thesemutvars),1)}];
end
for iallmutval = 1:cumprodnVals(end)
  for imutvar = 1:length(entry)
    if imutvar==1
      allmutvals{iallmutval}(imutvar) = entry{imutvar}{2}{1}(mod(iallmutval-1,nVals(imutvar))+1);
    else
      allmutvals{iallmutval}(imutvar) = entry{imutvar}{2}{1}(mod(floor((iallmutval-1)/cumprodnVals(imutvar-1)),nVals(imutvar))+1);          
    end
  end
end


str='';
for ivar=1:length(entry)
  underscoreind = strfind(entry{ivar}{1}{1},'_');
  if ~isempty(underscoreind)
    entry{ivar}{1}{1} = entry{ivar}{1}{1}(1:underscoreind-1);
  end
  if ~isempty(strfind(entry{ivar}{1}{1},'offm')) || ~isempty(strfind(entry{ivar}{1}{1},'offh')) || ~isempty(strfind(entry{ivar}{1}{1},'ehcn'))
    str = [str entry{ivar}{1}{1} ': ' repmat('+',allmutvals{myallmutval}(ivar)>=0,1) num2str(allmutvals{myallmutval}(ivar)) ' mV' coeff{1}];
  else
    str = [str entry{ivar}{1}{1} ': *' num2str(allmutvals{myallmutval}(ivar)) coeff{2}];
  end
  if ivar < length(entry)
    str = [str char(10)];
  end
end

      