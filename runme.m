% Runs the simulations needed for drawing Figure 1 and 2 and draws the figures.
% Assume scaling files are already calculated.
addpath('..');

cd kharche
runcontrol_kharche;
calcrates_kharche;
cd ..

cd severi
runcontrol_severi;
calcrates_severi;
cd ..

cd hay
% These commands may have to be run on command line in order to include the
% required paths (and parallelization is recommended anyway):
unix('nrnivmodl');
inds = [0 47 65 77 83 89 93 102 105];
for iind = 1:length(inds)
  unix(['python calcsteadystate.py ' num2str(inds(iind))]);      %each iind takes a couple of minutes to finish
  unix(['python findDCshortthreshold.py ' num2str(inds(iind))]); %each iind takes around 20 minutes to finish
  unix(['python calcifcurves.py ' num2str(inds(iind))]);         %each iind takes a couple of hours to finish
end
unix('python collectfig1.py');
unix('python collectfig2.py');
cd ..

cd almog
% These commands may have to be run on command line in order to include the
% required paths (and parallelization is recommended anyway):
unix('nrnivmodl');
inds = [0 47 65 77 83 89 93 102 105];
for iind = 1:length(inds)
  unix(['python calcsteadystate.py ' num2str(inds(iind))]);      %each iind takes around 40 minutes to finish
  unix(['python findDCshortthreshold.py ' num2str(inds(iind))]); %each iind takes around an hour to finish
  unix(['python calcifcurves.py ' num2str(inds(iind))]);         %each iind takes a couple of days to finish
end
unix('python collectfig1.py');
unix('python collectfig2.py');
cd ..

rmpath('..');

close all; drawfig1;
close all; drawfig2;
