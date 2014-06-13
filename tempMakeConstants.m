load smallCPDstimset.mat

numOrientations = 6;
dims = size(targetStim);

bandwidthList = [.5 1];

for iBandwidth = 1:length(bandwidthList)
    
    bandwidth = bandwidthList(iBandwidth);
    numLevels = maxLevel(dims,bandwidth);
    
    [freqRespsImag,freqRespsReal,pind]= ...
        makeQuadFRs(dims,numLevels,numOrientations,bandwidth);

    tempBlank = blankStim;
    stimBlank = tempBlank/mean2(tempBlank) - 1;
    [pyrBlank,pindBlank]=...
        buildQuadBands(stimBlank,freqRespsImag,freqRespsReal);
    energiesBlank = real(pyrBlank).^2 + imag(pyrBlank).^2;
    bandIndBlank = 1; % skip first and last images
    for iSF = 1:numLevels
        for iOR = 1:numOrientations
            bandIndBlank = bandIndBlank+1;
            EBlank{iBandwidth}(iSF,iOR,:,:) = reshape(energiesBlank(pyrBandIndices(pind,bandIndBlank)),dims(1),dims(2));
        end
    end
    smallFilterVals{iBandwidth}.freqRespsImag = freqRespsImag;
    smallFilterVals{iBandwidth}.freqRespsReal = freqRespsReal;
    smallFilterVals{iBandwidth}.pind = pind;
    smallFilterVals{iBandwidth}.EBlank = EBlank{iBandwidth};
    smallFilterVals{iBandwidth}.numLevels = numLevels;
    
    clear energiesBlank pyrBlank
end

save smallFilterVals.mat smallFilterVals

for iBandwidth = 1:length(bandwidthList)
    
    bandwidth = bandwidthList(iBandwidth);
    numLevels = maxLevel(dims,bandwidth);
    
    [freqRespsImag,freqRespsReal,pind]= ...
        makeQuadFRs(dims,numLevels,numOrientations,bandwidth);

    stimTarget = targetStim/mean2(targetStim) - 1;
    [pyrTarget,pindTarget]=...
        buildQuadBands(stimTarget,freqRespsImag,freqRespsReal);
    energiesTarget = real(pyrTarget).^2 + imag(pyrTarget).^2;
    bandIndTarget = 1; % skip first and last images
    for iSF = 1:numLevels
        for iOR = 1:numOrientations
            bandIndTarget = bandIndTarget+1;
            ETarget{iBandwidth}(iSF,iOR,:,:) = reshape(energiesTarget(pyrBandIndices(pind,bandIndTarget)),dims(1),dims(2));
        end
    end
    smallFilterVals{iBandwidth}.ETarget = ETarget{iBandwidth};
    
    clear energiesBlank pyrBlank
end

save smallFilterVals.mat smallFilterVals
