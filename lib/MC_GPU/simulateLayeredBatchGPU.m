function batch = simulateLayeredBatchGPU(nThis, GO, context)
%SIMULATELAYEREDBATCHGPU Advance one MC_GPU photon batch through the cloud.
%
% This GPU counterpart of the photon loop in MSMCSim keeps the validated
% legacy and packed state layouts as local implementations.

lidarCfg = context.lidarCfg;
cloudCfg = context.cloudCfg;
options = context.options;
stateLayout = context.stateLayout;
collectFullDiagnostics = context.collectFullDiagnostics;
maxScatterOrder = context.maxScatterOrder;
maxOrderEstimate = context.maxOrderEstimate;
precision = context.precision;

if strcmp(stateLayout, 'packed')
    batch = simulateLayeredBatchPacked(nThis, GO);
else
    batch = simulateLayeredBatchLegacy(nThis, GO);
end

    function batch = simulateLayeredBatchLegacy(nThis, GO)
        proto = GO.theta;
        rCl = lidarCfg.rBeam * 1e-3 + ...
            lidarCfg.divBeam * 1e-3 / 2 * cloudCfg.clBase;

        % Uniform launch over the complete circular beam cross-section.
        radius = cast(rCl, 'like', proto) .* sqrt(rand(nThis, 1, 'like', proto));
        azimuth = cast(2*pi, 'like', proto) .* rand(nThis, 1, 'like', proto);
        pPos = zeros(nThis, 3, 'like', proto);
        pPos(:, 1) = radius .* cos(azimuth);
        pPos(:, 2) = radius .* sin(azimuth);
        pPos(:, 3) = cast(cloudCfg.clBase, 'like', proto);
        pDir = pPos ./ sqrt(sum(pPos.^2, 2));

        SVIn = repmat(GO.SBeam, nThis, 1);
        SVIn = rotSVPhiGPU(SVIn, atan2(pDir(:, 2), pDir(:, 1)));

        pathLength = zeros(nThis, 1, 'like', proto);
        nScas = zeros(nThis, 1, 'int32', 'gpuArray');
        alive = true(nThis, 1, 'gpuArray');

        outSV = zeros(4, lidarCfg.nBins, 'like', proto);
        if collectFullDiagnostics
            outLayerSV = zeros( ...
                4, lidarCfg.nBins, GO.nLayers, 'like', proto);
            outByOrder = zeros( ...
                4, lidarCfg.nBins, maxOrderEstimate, 'like', proto);
            outOrderCount = zeros(maxOrderEstimate, 1, 'like', proto);
            outScatterCount = zeros(GO.nLayers, 1, 'like', proto);
        else
            outLayerSV = [];
            outByOrder = [];
            outOrderCount = [];
            outScatterCount = [];
        end

        for iStep = 1:maxScatterOrder
            if iStep == 1 || mod(iStep, 8) == 0
                if ~gather(any(alive))
                    break;
                end
            end

            activeGlobal = find(alive);
            if isempty(activeGlobal)
                break;
            end

            pos = pPos(activeGlobal, :);
            dir = pDir(activeGlobal, :);
            sv = SVIn(activeGlobal, :);
            len = pathLength(activeGlobal);
            ns = nScas(activeGlobal);

            [pos, isInMedium, ds, scatterLayer] = ...
                freeRunningGPU(pos, dir, GO, 'legacy');
            len = len + ds;
            if options.exactTimeWindowStop
                maxUsefulPath=cast(3e8*lidarCfg.binWidth*1e-9*lidarCfg.nBins,'like',proto);
                isInMedium=isInMedium & len<maxUsefulPath;
            end

            leaving = ~isInMedium;
            if any(leaving)
                alive(activeGlobal(leaving)) = false;
            end
            keep = find(isInMedium);
            if isempty(keep)
                pathLength(activeGlobal) = len;
                pPos(activeGlobal, :) = pos;
                continue;
            end

            pos = pos(keep, :);
            dir = dir(keep, :);
            sv = sv(keep, :);
            len = len(keep);
            ns = ns(keep);
            layer = scatterLayer(keep);
            globalIdx = activeGlobal(keep);

            survivesAbsorption = rand(numel(keep), 1, 'like', proto) <= ...
                GO.albedo(layer);
            killed = ~survivesAbsorption;
            if any(killed)
                alive(globalIdx(killed)) = false;
            end
            take = find(survivesAbsorption);
            if isempty(take)
                pPos(globalIdx, :) = pos;
                pDir(globalIdx, :) = dir;
                SVIn(globalIdx, :) = sv;
                pathLength(globalIdx) = len;
                continue;
            end

            pos = pos(take, :);
            dir = dir(take, :);
            sv = sv(take, :);
            len = len(take);
            ns = ns(take);
            layer = layer(take);
            globalIdx = globalIdx(take);

            if options.useRussianRoulette
                lowWeight = sv(:, 1) < options.energyThreshold;
                if any(lowWeight)
                    survivalProbability = sv(lowWeight, 1) ./ ...
                        cast(options.energyThreshold, 'like', proto);
                    lowLocal = find(lowWeight);
                    rouletteSurvives = rand(numel(lowLocal), 1, 'like', proto) <= ...
                        survivalProbability;
                    lowGlobal = globalIdx(lowWeight);
                    alive(lowGlobal(~rouletteSurvives)) = false;

                    surviveLocal = lowLocal(rouletteSurvives);
                    if ~isempty(surviveLocal)
                        survivingProbability = reshape(survivalProbability(rouletteSurvives), [], 1);
                        sv(surviveLocal, :) = sv(surviveLocal, :) ./ survivingProbability;
                    end
                    retain = ~lowWeight;
                    retain(surviveLocal) = true;

                    pos = pos(retain, :);
                    dir = dir(retain, :);
                    sv = sv(retain, :);
                    len = len(retain);
                    ns = ns(retain);
                    layer = layer(retain);
                    globalIdx = globalIdx(retain);
                end
            end

            if isempty(globalIdx)
                continue;
            end

            fovRadius = (lidarCfg.FOV * 1e-3 .* pos(:, 3) + ...
                lidarCfg.rTel * 1e-3 * 2) / 2;
            if strcmp(options.geometryMode,'mscart')
                fovRadius=tan(lidarCfg.FOV*1e-3/2).*pos(:,3);
            end
            inFOV = sqrt((pos(:, 1) - lidarCfg.distBeamTel * 1e-3).^2 + ...
                pos(:, 2).^2) < fovRadius;

            if any(inFOV)
                fovLocal = find(inFOV);
                fPos = pos(fovLocal, :);
                fDir = dir(fovLocal, :);
                fSV = sv(fovLocal, :);
                fLen = len(fovLocal);
                fNS = ns(fovLocal);
                fLayer = layer(fovLocal);
                if collectFullDiagnostics
                    outScatterCount = accumulateLayerCounts( ...
                        outScatterCount, fLayer, GO.nLayers, proto);
                end

                vectorToTelescope = [cast(lidarCfg.distBeamTel * 1e-3, 'like', proto) - fPos(:, 1), ...
                    -fPos(:, 2), -fPos(:, 3)];
                distanceSquared = sum(vectorToTelescope.^2, 2);
                pScaDir = vectorToTelescope ./ sqrt(distanceSquared);
                [scaAng, rotAng1, rotAng2] = scaRotAngGPU(fDir, pScaDir, 'legacy');
                iAng = nearestGridIndex(GO.theta, scaAng * 180 / pi);
                opticalIndex = fLayer + (iAng - 1) * GO.nLayers;

                layerBottom = GO.cloudBase + ...
                    (fLayer - 1) .* GO.layerThickness;
                opticalDepthCurrent = (fPos(:, 3) - layerBottom) .* ...
                    GO.muExt(fLayer);
                attenOD = GO.lowerTau(fLayer) + opticalDepthCurrent;
                if options.useLocalAtten
                    attenOD = opticalDepthCurrent;
                end
                probEnter = zeros(numel(fovLocal), 1, 'like', proto);
                receiverEligible = GO.theta(iAng) >= 90;
                if strcmp(options.forcedReceiverPathMode, 'allangles')
                    receiverEligible(:) = true;
                end
                if any(receiverEligible)
                    albedoBackward = GO.albedo(fLayer(receiverEligible));
                    s11Backward = GO.S11(opticalIndex(receiverEligible));
                    albedoBackward = albedoBackward(:);
                    s11Backward = s11Backward(:);
                    probEnter(receiverEligible) = albedoBackward .* ...
                        cast(pi * (lidarCfg.rTel * 1e-3)^2, 'like', proto) ./ ...
                        distanceSquared(receiverEligible) .* ...
                        exp(-attenOD(receiverEligible) ./ ...
                        max(abs(pScaDir(receiverEligible, 3)), eps(precision))) .* ...
                        s11Backward;
                end

                if strcmp(options.forcedReceiverRotation, 'vectorbasis')
                    [receiverInRotation, receiverOutRotation] = ...
                        recRotAng(fDir, pScaDir, proto);
                    forcedSV = applyMueller( ...
                        rotSVPhiGPU(fSV, receiverInRotation), ...
                        iAng, fLayer, GO);
                    forcedSV = rotSVPhiGPU(forcedSV, receiverOutRotation);
                else
                    receiverIncomingRotation = rotAng1;
                    if strcmp(options.forcedReceiverIncomingSign, 'donovan2015')
                        receiverIncomingRotation = -receiverIncomingRotation;
                    end
                    forcedSV = applyMueller( ...
                        rotSVPhiGPU(fSV, receiverIncomingRotation), ...
                        iAng, fLayer, GO);
                if strcmp(options.forcedReceiverRotation, 'donovan2015')
                    % Donovan et al. (2015), Eq. (A3): the outgoing-plane
                    % transform is L(pi-i2), which is equivalent to L(-i2)
                    % for Stokes Q/U.  scaRotAngGPU returns the supplementary
                    % angle (pi-i2), so the required argument is +rotAng2.
                    forcedSV = rotSVPhiGPU(forcedSV, rotAng2);
                else
                    forcedSV = rotSVPhiGPU(forcedSV, -rotAng2);
                end
                receiverAzimuth = atan2(pScaDir(:, 2), pScaDir(:, 1));
                if strcmp(options.forcedReceiverAzimuth, 'negative')
                    forcedSV = rotSVPhiGPU(forcedSV, -receiverAzimuth);
                elseif strcmp(options.forcedReceiverAzimuth, 'positive')
                    forcedSV = rotSVPhiGPU(forcedSV, receiverAzimuth);
                end
                end

                timeBin = (fLen + (fPos(:, 3) - cloudCfg.clBase) ./ ...
                    max(abs(pScaDir(:, 3)), eps(precision))) / 3e8 / ...
                    (lidarCfg.binWidth * 1e-9);
                if strcmp(options.geometryMode,'mscart')
                    timeBin=(fLen+sqrt(distanceSquared)-cloudCfg.clBase)/3e8/(lidarCfg.binWidth*1e-9);
                end
                iBin = ceil(timeBin);
                validBin = isfinite(iBin) & iBin >= 1 & iBin <= lidarCfg.nBins;
                if any(validBin)
                    pathFactor = 2.^(cast(fNS(validBin), 'like', proto) - 1);
                    if strcmp(options.forcedReceiverPathMode, 'allangles')
                        pathFactor(:) = cast(0.5,'like',proto); % Preserve the validated absolute signal convention
                    end
                    scale = probEnter(validBin) .* pathFactor;
                    finalSV = forcedSV(validBin, :) .* scale;
                    validIBin = iBin(validBin);
                    outSV = accumulateStokes(outSV, validIBin, finalSV, lidarCfg.nBins);
                    if collectFullDiagnostics
                        outLayerSV = accumulateLayerStokes( ...
                            outLayerSV, validIBin, fLayer(validBin), finalSV, ...
                            lidarCfg.nBins, GO.nLayers);

                        order = double(fNS(validBin)) + 1;
                        validOrder = order <= maxOrderEstimate;
                        if any(validOrder)
                            order = order(validOrder);
                            orderBin = validIBin(validOrder) + ...
                                (order - 1) * lidarCfg.nBins;
                            orderSV = finalSV(validOrder, :);
                            outByOrder = accumulateOrderStokes( ...
                                outByOrder, orderBin, orderSV, ...
                                lidarCfg.nBins, maxOrderEstimate);
                            outOrderCount = outOrderCount + accumarray( ...
                                order, ones(size(order), 'like', proto), ...
                                [maxOrderEstimate, 1], @sum, ...
                                cast(0, 'like', proto));
                        end
                    end
                end

                fSV = fSV .* (1 - probEnter);
                sv(fovLocal, :) = fSV;
            end

            rotAng1 = rand(numel(globalIdx), 1, 'like', proto) * (2 * pi);
            qRot = sv(:, 2) .* cos(2 * rotAng1) + ...
                sv(:, 3) .* sin(2 * rotAng1);
            iAng = sampleScatteringAngleGPU( ...
                GO.cdfS11, GO.cdfS12, sv(:, 1), qRot, ...
                layer, GO.nLayers, proto);
            scaAng = GO.theta(iAng) * (pi / 180);
            pScaDir = updateDirGPU(dir, rotAng1, scaAng, 'legacy');

            if strcmp(options.propagationRotation, 'vectorbasis')
                [propagationInRotation, propagationOutRotation] = ...
                    propagationRotationsVectorBasis(dir, pScaDir, proto);
                scatteredSV = applySampledMueller( ...
                    rotSVPhiGPU(sv, propagationInRotation), ...
                    iAng, layer, GO, ...
                    options.propagationWeightNormalization);
                sv = rotSVPhiGPU(scatteredSV, propagationOutRotation);
            else
                scatteredSV = applySampledMueller( ...
                    rotSVPhiGPU(sv, rotAng1), iAng, layer, GO, ...
                    options.propagationWeightNormalization);
                tempDen = sqrt(max(0, (1 - cos(scaAng).^2) .* ...
                    (1 - pScaDir(:, 3).^2)));
                meridianRotation = zeros(size(scaAng), 'like', proto);
                singular = abs(tempDen) < 1e-8;
                meridianRotation(singular) = pi / 2;
                regular = ~singular;
                if any(regular)
                    phiSign = 2 * cast(rotAng1(regular) > pi & ...
                        rotAng1(regular) < 2*pi, 'like', proto) - 1;
                    cosRotation = phiSign .* ...
                        (pScaDir(regular, 3) .* cos(scaAng(regular)) - ...
                        dir(regular, 3)) ./ tempDen(regular);
                    cosRotation = max(-1, min(1, cosRotation));
                    meridianRotation(regular) = acos(cosRotation);
                end
                sv = rotSVPhiGPU(scatteredSV, -meridianRotation);
            end
            ns = ns + 1;
            if collectFullDiagnostics
                outScatterCount = accumulateLayerCounts( ...
                    outScatterCount, layer, GO.nLayers, proto);
            end

            pPos(globalIdx, :) = pos;
            pDir(globalIdx, :) = pScaDir;
            SVIn(globalIdx, :) = sv;
            pathLength(globalIdx) = len;
            nScas(globalIdx) = ns;
        end

        batch = struct();
        batch.SVReturn = outSV;
        batch.truncatedPhotons = sum(alive);
        if collectFullDiagnostics
            batch.SVReturnLayer = outLayerSV;
            batch.SVByOrder = outByOrder;
            batch.orderPhotonCount = outOrderCount;
            batch.layerScatterCount = outScatterCount;
            batch.nScas = nScas;
        end
    end

    function batch = simulateLayeredBatchPacked(nThis, GO)
        proto = GO.theta;
        rCl = lidarCfg.rBeam * 1e-3 + ...
            lidarCfg.divBeam * 1e-3 / 2 * cloudCfg.clBase;

        % Uniform launch over the complete circular beam cross-section.
        radius = cast(rCl, 'like', proto) .* sqrt(rand(nThis, 1, 'like', proto));
        azimuth = cast(2*pi, 'like', proto) .* rand(nThis, 1, 'like', proto);
        pos = zeros(nThis, 3, 'like', proto);
        pos(:, 1) = radius .* cos(azimuth);
        pos(:, 2) = radius .* sin(azimuth);
        pos(:, 3) = cast(cloudCfg.clBase, 'like', proto);
        dir = pos ./ sqrt(sum(pos.^2, 2));

        sv = repmat(GO.SBeam, nThis, 1);
        sv = rotSVPhiGPU(sv, atan2(dir(:, 2), dir(:, 1)));
        len = zeros(nThis, 1, 'like', proto);
        ns = zeros(nThis, 1, 'int32', 'gpuArray');

        outSV = zeros(4, lidarCfg.nBins, 'like', proto);
        if collectFullDiagnostics
            outLayerSV = zeros( ...
                4, lidarCfg.nBins, GO.nLayers, 'like', proto);
            outByOrder = zeros( ...
                4, lidarCfg.nBins, maxOrderEstimate, 'like', proto);
            outOrderCount = zeros(maxOrderEstimate, 1, 'like', proto);
            outScatterCount = zeros(GO.nLayers, 1, 'like', proto);
            outNScasHist = zeros(maxScatterOrder + 1, 1, 'like', proto);
        else
            outLayerSV = [];
            outByOrder = [];
            outOrderCount = [];
            outScatterCount = [];
            outNScasHist = [];
        end

        for iStep = 1:maxScatterOrder
            if isempty(pos)
                break;
            end

            [pos, isInMedium, ds, scatterLayer] = ...
                freeRunningGPU(pos, dir, GO, 'packed');
            len = len + ds;
            if options.exactTimeWindowStop
                maxUsefulPath=cast(3e8*lidarCfg.binWidth*1e-9*lidarCfg.nBins,'like',proto);
                isInMedium=isInMedium & len<maxUsefulPath;
            end

            if collectFullDiagnostics
                outNScasHist = accumulateScatteringHistogram( ...
                    outNScasHist, ns(~isInMedium), maxScatterOrder, proto);
            end
            pos = pos(isInMedium, :);
            dir = dir(isInMedium, :);
            sv = sv(isInMedium, :);
            len = len(isInMedium);
            ns = ns(isInMedium);
            layer = scatterLayer(isInMedium);
            if isempty(pos)
                break;
            end

            survivesAbsorption = rand(size(pos, 1), 1, 'like', proto) <= ...
                GO.albedo(layer);
            if collectFullDiagnostics
                outNScasHist = accumulateScatteringHistogram( ...
                    outNScasHist, ns(~survivesAbsorption), ...
                    maxScatterOrder, proto);
            end
            pos = pos(survivesAbsorption, :);
            dir = dir(survivesAbsorption, :);
            sv = sv(survivesAbsorption, :);
            len = len(survivesAbsorption);
            ns = ns(survivesAbsorption);
            layer = layer(survivesAbsorption);
            if isempty(pos)
                break;
            end

            if options.useRussianRoulette
                lowWeight = sv(:, 1) < options.energyThreshold;
                lowLocal = reshape(find(lowWeight), [], 1);
                survivalProbability = sv(lowLocal, 1) ./ ...
                    cast(options.energyThreshold, 'like', proto);
                rouletteSurvives = rand(numel(lowLocal), 1, 'like', proto) <= ...
                    survivalProbability;
                killedLocal = lowLocal(~rouletteSurvives);
                if collectFullDiagnostics
                    outNScasHist = accumulateScatteringHistogram( ...
                        outNScasHist, ns(killedLocal), ...
                        maxScatterOrder, proto);
                end

                surviveLocal = lowLocal(rouletteSurvives);
                if ~isempty(surviveLocal)
                    survivingProbability = reshape(survivalProbability(rouletteSurvives), [], 1);
                    sv(surviveLocal, :) = sv(surviveLocal, :) ./ survivingProbability;
                end
                retain = ~lowWeight;
                retain(surviveLocal) = true;

                pos = pos(retain, :);
                dir = dir(retain, :);
                sv = sv(retain, :);
                len = len(retain);
                ns = ns(retain);
                layer = layer(retain);
                if isempty(pos)
                    break;
                end
            end

            fovRadius = (lidarCfg.FOV * 1e-3 .* pos(:, 3) + ...
                lidarCfg.rTel * 1e-3 * 2) / 2;
            if strcmp(options.geometryMode,'mscart')
                fovRadius=tan(lidarCfg.FOV*1e-3/2).*pos(:,3);
            end
            inFOV = sqrt((pos(:, 1) - lidarCfg.distBeamTel * 1e-3).^2 + ...
                pos(:, 2).^2) < fovRadius;

            fovLocal = reshape(find(inFOV), [], 1);
                fPos = pos(fovLocal, :);
                fDir = dir(fovLocal, :);
                fSV = sv(fovLocal, :);
                fLen = len(fovLocal);
                fNS = ns(fovLocal);
                fLayer = layer(fovLocal);
                if collectFullDiagnostics
                    outScatterCount = accumulateLayerCounts( ...
                        outScatterCount, fLayer, GO.nLayers, proto);
                end

                vectorToTelescope = [cast(lidarCfg.distBeamTel * 1e-3, 'like', proto) - fPos(:, 1), ...
                    -fPos(:, 2), -fPos(:, 3)];
                distanceSquared = sum(vectorToTelescope.^2, 2);
                pScaDir = vectorToTelescope ./ sqrt(distanceSquared);
                [scaAng, rotAng1, rotAng2] = ...
                    scaRotAngGPU(fDir, pScaDir, 'packed');
                iAng = nearestGridIndex(GO.theta, scaAng * 180 / pi);
                opticalIndex = fLayer + (iAng - 1) * GO.nLayers;

                layerBottom = GO.cloudBase + ...
                    (fLayer - 1) .* GO.layerThickness;
                opticalDepthCurrent = (fPos(:, 3) - layerBottom) .* ...
                    GO.muExt(fLayer);
                attenOD = GO.lowerTau(fLayer) + opticalDepthCurrent;
                if options.useLocalAtten
                    attenOD = opticalDepthCurrent;
                end
                probEnter = zeros(numel(fovLocal), 1, 'like', proto);
                receiverEligible = GO.theta(iAng) >= 90;
                if strcmp(options.forcedReceiverPathMode, 'allangles')
                    receiverEligible(:) = true;
                end
                backwardRows = reshape(find(receiverEligible), [], 1);
                albedoBackward = GO.albedo(fLayer(backwardRows));
                s11Backward = GO.S11(opticalIndex(backwardRows));
                albedoBackward = albedoBackward(:);
                s11Backward = s11Backward(:);
                probEnter(backwardRows) = albedoBackward .* ...
                    cast(pi * (lidarCfg.rTel * 1e-3)^2, 'like', proto) ./ ...
                    distanceSquared(backwardRows) .* ...
                    exp(-attenOD(backwardRows) ./ ...
                    max(abs(pScaDir(backwardRows, 3)), eps(precision))) .* ...
                    s11Backward;

                if strcmp(options.forcedReceiverRotation, 'vectorbasis')
                    [receiverInRotation, receiverOutRotation] = ...
                        recRotAng(fDir, pScaDir, proto);
                    forcedSV = applyMueller( ...
                        rotSVPhiGPU(fSV, receiverInRotation), ...
                        iAng, fLayer, GO);
                    forcedSV = rotSVPhiGPU(forcedSV, receiverOutRotation);
                else
                    receiverIncomingRotation = rotAng1;
                    if strcmp(options.forcedReceiverIncomingSign, 'donovan2015')
                        receiverIncomingRotation = -receiverIncomingRotation;
                    end
                    forcedSV = applyMueller( ...
                        rotSVPhiGPU(fSV, receiverIncomingRotation), ...
                        iAng, fLayer, GO);
                    if strcmp(options.forcedReceiverRotation, 'donovan2015')
                        forcedSV = rotSVPhiGPU(forcedSV, rotAng2);
                    else
                        forcedSV = rotSVPhiGPU(forcedSV, -rotAng2);
                    end
                    receiverAzimuth = atan2(pScaDir(:, 2), pScaDir(:, 1));
                    if strcmp(options.forcedReceiverAzimuth, 'negative')
                        forcedSV = rotSVPhiGPU(forcedSV, -receiverAzimuth);
                    elseif strcmp(options.forcedReceiverAzimuth, 'positive')
                        forcedSV = rotSVPhiGPU(forcedSV, receiverAzimuth);
                    end
                end

                timeBin = (fLen + (fPos(:, 3) - cloudCfg.clBase) ./ ...
                    max(abs(pScaDir(:, 3)), eps(precision))) / 3e8 / ...
                    (lidarCfg.binWidth * 1e-9);
                if strcmp(options.geometryMode,'mscart')
                    timeBin=(fLen+sqrt(distanceSquared)-cloudCfg.clBase)/3e8/(lidarCfg.binWidth*1e-9);
                end
                iBin = ceil(timeBin);
                validBin = isfinite(iBin) & iBin >= 1 & iBin <= lidarCfg.nBins;
                validRows = reshape(find(validBin), [], 1);
                pathFactor = 2.^(cast(fNS(validRows), 'like', proto) - 1);
                if strcmp(options.forcedReceiverPathMode, 'allangles')
                    pathFactor(:) = cast(0.5,'like',proto); % Preserve the validated absolute signal convention
                end
                scale = probEnter(validRows) .* pathFactor;
                scale = scale(:);
                finalSV = forcedSV(validRows, :) .* scale;
                validIBin = iBin(validRows);
                validIBin = validIBin(:);
                outSV = accumulateStokes(outSV, validIBin, finalSV, lidarCfg.nBins);
                if collectFullDiagnostics
                    outLayerSV = accumulateLayerStokes( ...
                        outLayerSV, validIBin, fLayer(validRows), finalSV, ...
                        lidarCfg.nBins, GO.nLayers);

                    order = double(fNS(validRows)) + 1;
                    validOrder = order <= maxOrderEstimate;
                    order = order(validOrder);
                    orderBin = validIBin(validOrder) + ...
                        (order - 1) * lidarCfg.nBins;
                    orderSV = finalSV(validOrder, :);
                    outByOrder = accumulateOrderStokes( ...
                        outByOrder, orderBin, orderSV, ...
                        lidarCfg.nBins, maxOrderEstimate);
                    outOrderCount = outOrderCount + accumarray( ...
                        order, ones(size(order), 'like', proto), ...
                        [maxOrderEstimate, 1], @sum, ...
                        cast(0, 'like', proto));
                end

                fSV = fSV .* (1 - probEnter);
                sv(fovLocal, :) = fSV;

            rotAng1 = rand(size(pos, 1), 1, 'like', proto) * (2 * pi);
            qRot = sv(:, 2) .* cos(2 * rotAng1) + ...
                sv(:, 3) .* sin(2 * rotAng1);
            iAng = sampleScatteringAngleGPU( ...
                GO.cdfS11, GO.cdfS12, sv(:, 1), qRot, ...
                layer, GO.nLayers, proto);
            scaAng = GO.theta(iAng) * (pi / 180);
            pScaDir = updateDirGPU(dir, rotAng1, scaAng, 'packed');

            if strcmp(options.propagationRotation, 'vectorbasis')
                [propagationInRotation, propagationOutRotation] = ...
                    propagationRotationsVectorBasis(dir, pScaDir, proto);
                scatteredSV = applySampledMueller( ...
                    rotSVPhiGPU(sv, propagationInRotation), ...
                    iAng, layer, GO, ...
                    options.propagationWeightNormalization);
                sv = rotSVPhiGPU(scatteredSV, propagationOutRotation);
            else
                scatteredSV = applySampledMueller( ...
                    rotSVPhiGPU(sv, rotAng1), iAng, layer, GO, ...
                    options.propagationWeightNormalization);
                tempDen = sqrt(max(0, (1 - cos(scaAng).^2) .* ...
                    (1 - pScaDir(:, 3).^2)));
                meridianRotation = zeros(size(scaAng), 'like', proto);
                singular = abs(tempDen) < 1e-8;
                meridianRotation(singular) = pi / 2;
                regular = ~singular;
                regularRows = reshape(find(regular), [], 1);
                phiSign = 2 * cast(rotAng1(regularRows) > pi & ...
                    rotAng1(regularRows) < 2*pi, 'like', proto) - 1;
                cosRotation = phiSign .* ...
                    (pScaDir(regularRows, 3) .* cos(scaAng(regularRows)) - ...
                    dir(regularRows, 3)) ./ tempDen(regularRows);
                cosRotation = max(-1, min(1, cosRotation));
                meridianRotation(regularRows) = acos(cosRotation);
                sv = rotSVPhiGPU(scatteredSV, -meridianRotation);
            end
            dir = pScaDir;
            ns = ns + 1;
            if collectFullDiagnostics
                outScatterCount = accumulateLayerCounts( ...
                    outScatterCount, layer, GO.nLayers, proto);
            end
        end

        truncatedCount = size(pos, 1);
        if collectFullDiagnostics
            outNScasHist = accumulateScatteringHistogram( ...
                outNScasHist, ns, maxScatterOrder, proto);
        end

        batch = struct();
        batch.SVReturn = outSV;
        batch.truncatedPhotons = cast(truncatedCount, 'like', proto);
        if collectFullDiagnostics
            batch.SVReturnLayer = outLayerSV;
            batch.SVByOrder = outByOrder;
            batch.orderPhotonCount = outOrderCount;
            batch.layerScatterCount = outScatterCount;
            batch.nScasHist = outNScasHist;
        end
    end

end

function out = applyMueller(in, idx, layer, G)
opticalIndex = layer + (idx - 1) * G.nLayers;
s11 = G.S11(opticalIndex);
s12 = G.S12(opticalIndex);
s33 = G.S33(opticalIndex);
s34 = G.S34(opticalIndex);
s11 = s11(:);
s12 = s12(:);
s33 = s33(:);
s34 = s34(:);
out = in;
out(:, 1) = (s11 .* in(:, 1) + s12 .* in(:, 2)) ./ s11;
out(:, 2) = (s12 .* in(:, 1) + s11 .* in(:, 2)) ./ s11;
out(:, 3) = (s33 .* in(:, 3) + s34 .* in(:, 4)) ./ s11;
out(:, 4) = (-s34 .* in(:, 3) + s33 .* in(:, 4)) ./ s11;
end

function out = applySampledMueller(in, idx, layer, G, normalizationMode)
%APPLYSAMPLEDMUELLER Account for the proposal used to sample direction.
% The legacy implementation divides by S11.  When the proposal density is
% proportional to S11*I + S12*Q, importance sampling requires division by
% that same polarized density, including its discrete CDF normalization.
opticalIndex = layer + (idx - 1) * G.nLayers;
s11 = G.S11(opticalIndex);
s12 = G.S12(opticalIndex);
s33 = G.S33(opticalIndex);
s34 = G.S34(opticalIndex);
s11 = s11(:);
s12 = s12(:);
s33 = s33(:);
s34 = s34(:);

denominator = s11;
if strcmp(normalizationMode, 'samplingpdf')
    nTheta = size(G.cdfS11, 2);
    lastIndex = layer + (nTheta - 1) * G.nLayers;
    cdf11End = G.cdfS11(lastIndex);
    cdf12End = G.cdfS12(lastIndex);
    cdf11End = cdf11End(:);
    cdf12End = cdf12End(:);
    proposalDensity = s11 .* in(:, 1) + s12 .* in(:, 2);
    proposalTotal = in(:, 1) .* cdf11End + in(:, 2) .* cdf12End;
    denominator = proposalDensity .* cdf11End ./ ...
        max(proposalTotal, eps(classUnderlying(proposalTotal)));
    denominator = max(denominator, eps(classUnderlying(denominator)));
end

out = in;
out(:, 1) = (s11 .* in(:, 1) + s12 .* in(:, 2)) ./ denominator;
out(:, 2) = (s12 .* in(:, 1) + s11 .* in(:, 2)) ./ denominator;
out(:, 3) = (s33 .* in(:, 3) + s34 .* in(:, 4)) ./ denominator;
out(:, 4) = (-s34 .* in(:, 3) + s33 .* in(:, 4)) ./ denominator;
if strcmp(normalizationMode,'p11')
    correction=G.proposalCorrection(opticalIndex);
    out=out.*correction(:);
end
end

function [rotationIn, rotationOut] = ...
    recRotAng(directionIn, directionOut, proto)
%RECROTANG Build the two Stokes rotations from vectors.
%
% The incident packet Stokes vector is expressed in its meridian-plane
% basis.  The returned signal is expressed in the receiver's fixed x-axis
% basis.  Constructing both rotations from projected unit vectors avoids
% ambiguity from supplementary-angle formulae near lidar backscatter.

n = cross(directionIn, directionOut, 2);
nNorm = sqrt(sum(n.^2, 2));
regular = nNorm > cast(1e-10, 'like', proto);
n = n ./ max(nNorm, cast(1e-10, 'like', proto));

eMeridian = meridianParallelAxis(directionIn, proto);
eDetector = fixedXAxisProjection(directionOut, proto);

pIn = cross(n, directionIn, 2);
pOut = cross(n, directionOut, 2);
pIn = normalizeRows(pIn, proto);
pOut = normalizeRows(pOut, proto);

% In the collinear limit the scattering plane is arbitrary.  Use the
% incident meridian axis as that arbitrary plane; this preserves the
% expected zero depolarization for single 180-degree scattering by spheres.
if any(~regular)
    pIn(~regular, :) = eMeridian(~regular, :);
    pOut(~regular, :) = eMeridian(~regular, :);
end

rotationIn = signedAxisAngle(eMeridian, pIn, directionIn);
rotationOut = signedAxisAngle(pOut, eDetector, directionOut);
end

function [rotationIn, rotationOut] = ...
    propagationRotationsVectorBasis(directionIn, directionOut, proto)
%PROPAGATIONROTATIONSVECTORBASIS Rotate between meridian/scattering planes.
n = cross(directionIn, directionOut, 2);
nNorm = sqrt(sum(n.^2, 2));
regular = nNorm > cast(1e-10, 'like', proto);
n = n ./ max(nNorm, cast(1e-10, 'like', proto));

eMeridianIn = meridianParallelAxis(directionIn, proto);
eMeridianOut = meridianParallelAxis(directionOut, proto);
pIn = normalizeRows(cross(n, directionIn, 2), proto);
pOut = normalizeRows(cross(n, directionOut, 2), proto);
if any(~regular)
    pIn(~regular, :) = eMeridianIn(~regular, :);
    pOut(~regular, :) = eMeridianOut(~regular, :);
end
rotationIn = signedAxisAngle(eMeridianIn, pIn, directionIn);
rotationOut = signedAxisAngle(pOut, eMeridianOut, directionOut);
end

function axis = meridianParallelAxis(direction, proto)
horizontal = sqrt(direction(:, 1).^2 + direction(:, 2).^2);
axis = zeros(size(direction), 'like', proto);
regular = horizontal > cast(1e-10, 'like', proto);
axis(regular, 1) = direction(regular, 3) .* ...
    direction(regular, 1) ./ horizontal(regular);
axis(regular, 2) = direction(regular, 3) .* ...
    direction(regular, 2) ./ horizontal(regular);
axis(regular, 3) = -horizontal(regular);
axis(~regular, 1) = 1;
end

function axis = fixedXAxisProjection(direction, proto)
axis = zeros(size(direction), 'like', proto);
axis(:, 1) = 1;
axis = axis - direction .* direction(:, 1);
axis = normalizeRows(axis, proto);
end

function rows = normalizeRows(rows, proto)
rowNorm = sqrt(sum(rows.^2, 2));
rows = rows ./ max(rowNorm, cast(1e-12, 'like', proto));
end

function angle = signedAxisAngle(axisFrom, axisTo, propagationDirection)
crossTerm = cross(axisFrom, axisTo, 2);
sineTerm = sum(propagationDirection .* crossTerm, 2);
cosineTerm = sum(axisFrom .* axisTo, 2);
angle = atan2(sineTerm, cosineTerm);
end

function idx = nearestGridIndex(grid, x)
lo = ones(size(x), 'like', x);
hi = numel(grid) * ones(size(x), 'like', x);
for k = 1:ceil(log2(double(numel(grid))))
    mid = floor((lo + hi) / 2);
    moveLeft = grid(mid) >= x;
    hi(moveLeft) = mid(moveLeft);
    lo(~moveLeft) = mid(~moveLeft) + 1;
end
idx = min(lo, numel(grid));
previous = max(1, idx - 1);
usePrevious = abs(grid(previous) - x) < abs(grid(idx) - x);
idx(usePrevious) = previous(usePrevious);
end

function out = accumulateStokes(out, bins, values, nBins)
zero = cast(0, 'like', values);
for iComponent = 1:4
    out(iComponent, :) = out(iComponent, :) + ...
        accumarray(bins, values(:, iComponent), [nBins, 1], @sum, zero).';
end
end

function out = accumulateOrderStokes(out, orderBins, values, nBins, nOrders)
zero = cast(0, 'like', values);
for iComponent = 1:4
    accumulated = accumarray(orderBins, values(:, iComponent), ...
        [nBins * nOrders, 1], @sum, zero);
    out(iComponent, :, :) = out(iComponent, :, :) + ...
        reshape(accumulated, [1, nBins, nOrders]);
end
end

function out = accumulateLayerStokes( ...
    out, bins, layers, values, nBins, nLayers)
zero = cast(0, 'like', values);
layerBins = bins + (layers - 1) * nBins;
for iComponent = 1:4
    accumulated = accumarray(layerBins, values(:, iComponent), ...
        [nBins * nLayers, 1], @sum, zero);
    out(iComponent, :, :) = out(iComponent, :, :) + ...
        reshape(accumulated, [1, nBins, nLayers]);
end
end

function out = accumulateLayerCounts(out, layers, nLayers, proto)
out = out + accumarray(layers, ones(size(layers), 'like', proto), ...
    [nLayers, 1], @sum, cast(0, 'like', proto));
end

function out = accumulateScatteringHistogram(out, nScas, maxOrder, proto)
if isempty(nScas)
    return;
end
bins = min(double(nScas), maxOrder) + 1;
out = out + accumarray(bins, ones(size(bins), 'like', proto), ...
    [maxOrder + 1, 1], @sum, cast(0, 'like', proto));
end
