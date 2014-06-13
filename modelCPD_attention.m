function results = modelCPD_attention(indexVal,inVals)

if nargin==0
    error('modelCPD_4D_HPC: need input argument')
end

keyboard
numOrientations = 6;
numSimulatedSystems = 10;

systemList = 1:numSimulatedSystems;
attentionFieldSizeList = [0 2 4 16]
bandwidthList = [.5 1];
sigmaList = linspace(1e-6,1e-4,3);
numeratorWidthList = [.5 1 2];
denominatorFactorList = [2 4 8];

parameterSpace = inVals.parameterSpace;
numParams = length(parameterSpace);
currentParams = parameterSpace(indexVal,:);

curAxWidth = currentParams(1);
curBandwidth = currentParams(2);
curSigma = currentParams(3);
curNumWidth = currentParams(4);
curDenomFactor = currentParams(5);

blankStim = inVals.stimSet.blankStim;
targetStim = inVals.stimSet.targetStim;

dims = size(blankStim);

switch curBandwidth
    case .5
        freqRespsImag = inVals.filtersVals{1}.freqRespsImag;
        freqRespsReal = inVals.filtersVals{1}.freqRespsReal;
        pind = inVals.filtersVals{1}.pind;
        blank = inVals.filtersVals{1}.EBlank;
        target = inVals.filtersVals{1}.ETarget;
        numLevels = inVals.filtersVals{1}.numLevels;

    case 1
        freqRespsImag = inVals.filtersVals{2}.freqRespsImag;
        freqRespsReal = inVals.filtersVals{2}.freqRespsReal;
        pind = inVals.filtersVals{2}.pind;
        blank = inVals.filtersVals{2}.EBlank;
        target = inVals.filtersVals{2}.ETarget;
        numLevels = inVals.filtersVals{2}.numLevels;

end

clear inVals

if(~curAxWidth) % if AxWidth = 0, no attention
    Eb = blank.^2;
    Et = target.^2;
else % multiply by attention field and then filter
     attnGainX = mkGaussian(dims,curAxWidth);
     attnField = repmat(attnGainX,[1,1,numLevels,numOrientations]);
     attnField = shiftdim(attnField,2);
     Eb = (blank.*attnField).^2;
     Et = (target.*attnField).^2;
end

% -- *** -- The subbands need blur size corresponding to their SF --***--
blurValTemp = 4*[1, 2, 4, 8, 16, 32, 64, linspace(100,500,max(1,numLevels-7))];
blurVal = blurValTemp(1:numLevels);

boundProjection = NaN*ones(numLevels); 

% normalize across all orientations (start with mean)
% first, blank:
I_OR_blank = mean(Eb,2);
IBlank = repmat(I_OR_blank,[1,numOrientations,1,1]);
EBlank = Eb./(IBlank+curSigma);

I_xy_blank = NaN*ones(numLevels, numOrientations, dims(1),dims(2));

% normalize across 2D space
for iband = 1:numLevels
    
    ker = mkGaussian([201,1],blurVal(iband));
    
    for iOR = 1:numOrientations
        
        EBlank_temp = squeeze(EBlank(iband,iOR,:,:));
        temp = conv2(ker,ker,EBlank_temp,'same');
        
        I_xy_blank(iband,iOR,:,:) = temp;
        temp = [];
        EBlank_temp = [];
        
    end
end

EBlank = EBlank./(I_xy_blank+curSigma);

% then, target:
I_OR_target = mean(Et,2);
ITarget = repmat(I_OR_target,[1,numOrientations,1,1]);
ETarget = Et./(ITarget+curSigma);

I_xy_target = NaN*ones(numLevels, numOrientations, dims(1),dims(2));

% normalize across 2D space
for iband = 1:numLevels
    
    ker = mkGaussian([201,1],blurVal(iband));
    
    for iOR = 1:numOrientations
        
        ETarget_temp = squeeze(ETarget(iband,iOR,:,:));
        temp = conv2(ker,ker,ETarget_temp,'same');
        
        I_xy_target(iband,iOR,:,:) = temp;
        temp = [];
        ETarget_temp = [];
        
    end
end

ETarget = ETarget./(I_xy_target+curSigma);



% ------------------------------------------------------------------------
% (4) model SF normalization
% ------------------------------------------------------------------------


% we have numSys (10) systems, which will have weighted sums of the SF
% level bands available to them, progressively leaving out high SF bands as
% move to periphery.

% Wbase is 10xnumLevels. The first row is a system with equal input
% from all SF channels. The last row is a system missing
% high frequency channels.
% The width of the convolution function determines how the channels are
% pooled

% The inhibitory denominator pools from more neurons.


firstValList = linspace(0,9,numSimulatedSystems);

for iSys = 1:numSimulatedSystems

    baseW = linspace(10.^(-firstValList(iSys)),1,numLevels);
    
    poolFilterNum = pdf('normal',linspace(1,numLevels,numLevels/2),iSys,curNumWidth);
    temp_num = conv(baseW,poolFilterNum);
    sumVal = sum(temp_num(1:numLevels));
    W_numerator(iSys,:) = temp_num(1:numLevels)./sumVal;
    clear poolFilterNum temp_num sumVal

    poolFilterDenom = pdf('normal',linspace(1,numLevels,numLevels/2),iSys,curNumWidth*curDenomFactor);
    temp_denom = conv(baseW,poolFilterDenom);
    sumVal = sum(temp_denom(1:numLevels));
    W_denominator(iSys,:) = temp_denom(1:numLevels)./sumVal;
    clear poolFilterDenom temp_denom sumVal
end


for iSystem = 1:numSimulatedSystems
    % model systems by subband weighting
    Eband_weighting = repmat(transpose(W_numerator(iSystem,:)),[1,numOrientations, dims(1),dims(2)]);
    Iband_weighting = repmat(transpose(W_numerator(iSystem,:)),[1,numOrientations, dims(1),dims(2)]);
    
    ETarget_Sys{iSystem} = ETarget.*Eband_weighting;
    ITarget_Sys{iSystem} = ETarget.*Iband_weighting;
    RTarget_Sys{iSystem} = ETarget_Sys{iSystem}./(ITarget_Sys{iSystem}+curSigma);
    
    % ------------------------------------------------------------------------
    % (5)   Create similar matrices for blank stim, CW and CCW
   
    
    
    EBlank_Sys{iSystem} = EBlank.*Eband_weighting;
    IBlank_Sys{iSystem} = EBlank.*Iband_weighting;
    RBlank_Sys{iSystem} = EBlank_Sys{iSystem}./(IBlank_Sys{iSystem} + curSigma);    
    
    boundary_Sys = RTarget_Sys{iSystem}(:)-RBlank_Sys{iSystem}(:);
    boundary_Sys = boundary_Sys/norm(boundary_Sys);
    
    % ------------------------------------------------------------------------
    % (6) estimate d' by projecting against the difference matrix
    
    targetProjection{iSystem} = boundary_Sys(:)'*RTarget_Sys{iSystem}(:);
    blankProjection{iSystem} = boundary_Sys(:)'*RBlank_Sys{iSystem}(:);
    
    dMeasure{iSystem} = abs(targetProjection{iSystem}-blankProjection{iSystem});
 
end % going through systems

results.targetProjection = targetProjection;
results.blankProjection = blankProjection;
results.dMeasure = dMeasure;
results.parameterSetIndexVal = indexVal;
results.parameterSet = parameterSpace(indexVal,:);

% -------------------------------
%-----------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

% SUBFUNCTIONS


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function boundary = getBoundary(A,B)

boundary = A-B;
boundary = boundary/norm(boundary(:));


function originalNormalizationCode

% Stimuulus drive
Eraw = conv2sepYcirc(stimulus,ExKernel,EthetaKernel) + baselineMod;
Emax = max(Eraw(:));
E = attnGain .* Eraw;

% Suppressive drive
I = conv2sepYcirc(E,IxKernel,IthetaKernel);
Imax = max(I(:));

% Normalization
R = E ./ (I + sigma) + baselineUnmod;



