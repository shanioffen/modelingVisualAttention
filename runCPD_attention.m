function runCPD_attention

matPath = '';

addpath('~/Dropbox/NYU/matlab/pyrTools/')
addpath('~/Dropbox/NYU/matlab/attentionModel/steer')

numOrientations = 6;
numSimulatedSystems = 10;

systemList = 1:numSimulatedSystems;
attentionFieldSizeList = [0 2 8 32 128];
bandwidthList = [.5 1];
sigmaList = linspace(1e-6,1e-4,3);
numeratorWidthList = [.5 1 2];
denominatorFactorList = [2 4 8];

parameterSpace = ...
  allcomb(attentionFieldSizeList, bandwidthList,sigmaList,numeratorWidthList,denominatorFactorList);

numParams = length(parameterSpace);
load([matPath 'smallCPDstimset.mat']);

load([matPath 'smallFilterVals.mat']);

inVals.filtersVals = smallFilterVals;

inVals.stimSet.blankStim = blankStim;
inVals.stimSet.targetStim = targetStim;
inVals.parameterSpace = parameterSpace;

clear smallFilterVals

results = [];

for iRun = 1%:numParams
    results{iRun} = modelCPD_attention(iRun,inVals);
    save('-v7.3',['intermediateResults' num2str(iRun) '.mat'],'results')
end

save -v7.3 allResultsCPD_updated.mat
