nEnsemble = 600;

oceanInitialObs = 31;
windInitialObs = 31;

load('approx_params600');
load("windou.mat")

dto = 1/108;

dtf = 30;

stepsPerIt = 3600*24/itsPerDay/dtf;
floePerOcn = round(dto*3600*24/dtf);
windStep = round(1/216*3600*24/dtf);

Xo = ocean.Xo;
Yo = ocean.Yo;

xw = linspace(xMin, xMax, 128);
yw = linspace(yMin, yMax, 128);
[Xw, Yw] = meshgrid(xw, yw); 

sigobs = 250;
sigalpha = pi/36;
sigvel = 250;

floeEnsemble = repmat(transpose(FloeCopy), 1, nEnsemble);
xEnsemble = nan(nFloes, nEnsemble, nIterations);
yEnsemble = nan(nFloes, nEnsemble, nIterations);
alphaEnsemble = nan(nFloes, nEnsemble, nIterations);

xArray = nan(nFloes, nEnsemble);
yArray = nan(nFloes, nEnsemble);

alphaArray = nan(nFloes, nEnsemble);
ksiArray = nan(nFloes, nEnsemble);
dxArray = nan(nFloes, nEnsemble);
dyArray = nan(nFloes, nEnsemble);
dalphaArray = nan(nFloes, nEnsemble);
uArray = nan(nFloes, nEnsemble);
vArray = nan(nFloes, nEnsemble);
duArray = nan(nFloes, nEnsemble);
dvArray = nan(nFloes, nEnsemble);
dksiArray = nan(nFloes, nEnsemble);

if ~exist('thicknessArray', 'var')
    thicknessArray = gamrnd(4, 0.5, nFloes, nEnsemble)+0.5;
    for iFloe = 1:nFloes
        for iEnsemble = 1:nEnsemble
            thickness = thicknessArray(iFloe, iEnsemble);
            floeEnsemble(iFloe, iEnsemble).mass = floeEnsemble(iFloe, iEnsemble).mass/floeEnsemble(iFloe, iEnsemble).thickness*thickness;
            floeEnsemble(iFloe, iEnsemble).inertia_moment = floeEnsemble(iFloe, iEnsemble).inertia_moment/floeEnsemble(iFloe, iEnsemble).thickness*thickness;
            floeEnsemble(iFloe, iEnsemble).thickness = thickness;
            
            floeEnsemble(iFloe, iEnsemble).Ui = randn(1)*1;
            floeEnsemble(iFloe, iEnsemble).Vi = randn(1)*1;

            floeEnsemble(iFloe, iEnsemble).ksi_ice = randn(1)*2e-6;
        end
    end
end



maxMode = 11;
nEnsModes = length(-maxMode:maxMode);
modeIdx = [1:maxMode+1, nModes-maxMode+1:nModes];


phArray = nan(nEnsModes, nEnsModes, nEnsemble);
for iEnsemble = 1:nEnsemble
    ph = eqMean + sqrt(eqVar).*(randn(nModes, nModes) + 1i*randn(nModes, nModes))/sqrt(2);
    phArray(:, :, iEnsemble) = ph(modeIdx, modeIdx);
end
phInitial = phArray;

maxWindMode = 4;
nWindModes = length(-maxWindMode:maxWindMode);
windModeIdx = [1:maxWindMode+1, nWindModes-maxWindMode+1:nWindModes];

windHatArray = nan(nWindModes, nWindModes, nEnsemble);
for iEnsemble = 1:nEnsemble
    wind = windEqMean + sqrt(windEqVar).*(randn(128, 128) + 1i*randn(128, 128))/sqrt(2);
    windHatArray(:, :, iEnsemble) = wind(windModeIdx, windModeIdx);
end
windHatInitial = windHatArray;

phEnsMean = nan(nEnsModes, nEnsModes, nObservations+1);
phEnsMean(:, :, 1) = mean(phArray, 3);

thicknessEnsMean = nan(nFloes, nObservations+1);
thicknessEnsStd = nan(nFloes, nObservations+1);

thicknessEnsMean(:, 1) = mean(thicknessArray, 2);
thicknessEnsStd(:, 1) = std(thicknessArray, 0, 2);

fprintf('Obs:')
for iObservation = 1:nObservations-1

    thicknessEnsMean(:, iObservation+1) = mean(thicknessArray, 2);
    thicknessEnsStd(:, iObservation+1) = std(thicknessArray, 0, 2);

    phEnsMean(:, :, iObservation+1) = mean(phArray, 3);

    iterationIdx = obsIdx(iObservation)+1:obsIdx(iObservation+1); % Needed for parfor

    activeFloes = obsFloesCopy(:, iObservation) & obsFloesCopy(:, iObservation+1);
    if nnz(activeFloes) == 0
        fprintf(' %d', iObservation);
        continue
    end

    parfor iEnsemble = 1:nEnsemble
        
        for iIteration = iterationIdx
    
            windCoeff = 1;



            ph = zeros(nModes);
            ph(modeIdx, modeIdx) = squeeze(phArray(:, :, iEnsemble));


            windHat = zeros(128);
            windHat(windModeIdx, windModeIdx) = squeeze(windHatArray(:, :, iEnsemble));
            
            for iStep = 1:stepsPerIt
                
                if mod(iStep+iIteration, floePerOcn) == 0
                    dW = sqrt(dto)*(randn(nModes) + 1i*randn(nModes))/sqrt(2);
                    dph = ((-abs(damping) + phase*1i).*ph + forcing)*dto + noise.*dW;
                    ph = ph + dph;                   
                end
                
                dtWind = dto*24;
                if mod(iStep+iIteration, floePerOcn) == 0
                    dWWind = sqrt(dtWind)*(randn(128) + 1i*randn(128))/sqrt(2);

                    dWind = ((-abs(windDamping) + windPhase*1i).*windHat + windForcing)*dtWind + windNoise.*dWWind;

                    windHat = windHat + dWind;                   
                end

                [u,v] = caluv(ph,k,l,trunc); 
                
                oceanTmp = struct('Uocn', u*1e3/0.86e5, 'Vocn', v*1e3/0.86e5, 'Xo', Xo, 'Yo', Yo);

                wind = ifft2(windHat);
                winds = struct('Xw', Xw, 'Yw', Yw, 'Uw', windCoeff*real(wind), 'Vw', windCoeff*imag(wind));
                
                floeEnsemble(:, iEnsemble) = calc_trajectory_multi_p(dtf,oceanTmp, winds, floeEnsemble(:, iEnsemble), floeShapes, activeFloes);
            end
            
            phArray(:, :, iEnsemble) = ph(modeIdx, modeIdx);

            windHatArray(:, :, iEnsemble) = windHat(windModeIdx, windModeIdx);
            
            xEnsemble(:, iEnsemble, iIteration) = [floeEnsemble(:, iEnsemble).Xi].';
            yEnsemble(:, iEnsemble, iIteration) = [floeEnsemble(:, iEnsemble).Yi].';
            alphaEnsemble(:, iEnsemble, iIteration) = [floeEnsemble(:, iEnsemble).alpha_i].';
            
            xArray(:, iEnsemble) = [floeEnsemble(:, iEnsemble).Xi].';
            yArray(:, iEnsemble) = [floeEnsemble(:, iEnsemble).Yi].';
            
            alphaArray(:, iEnsemble) = [floeEnsemble(:, iEnsemble).alpha_i].';
            ksiArray(:, iEnsemble) = [floeEnsemble(:, iEnsemble).ksi_ice].';
            uArray(:, iEnsemble) = [floeEnsemble(:, iEnsemble).Ui].';
            vArray(:, iEnsemble) = [floeEnsemble(:, iEnsemble).Vi].';
            dxArray(:, iEnsemble) = [floeEnsemble(:, iEnsemble).dXi_p].';
            dyArray(:, iEnsemble) = [floeEnsemble(:, iEnsemble).dYi_p].';
            dalphaArray(:, iEnsemble) = [floeEnsemble(:, iEnsemble).dalpha_i_p].';
            duArray(:, iEnsemble) = [floeEnsemble(:, iEnsemble).dUi_p].';
            dvArray(:, iEnsemble) = [floeEnsemble(:, iEnsemble).dVi_p].';
            dksiArray(:, iEnsemble) = [floeEnsemble(:, iEnsemble).dksi_ice_p].';
            
        end
    end

    if iObservation <= oceanInitialObs
        phInitial = phArray;
    end

    if iObservation <= windInitialObs
        windHatInitial = windHatArray;
    end
    

    nObsFloes = nnz(observedFloes(:, iObservation+1));
    obsFloesIdx = find(observedFloes(:, iObservation+1));



    
    if nObsFloes == 0
        fprintf(' %d', iObservation);
        continue
    end
    
    
    nObsVar = 3;
    dObserved = nObsVar*nObsFloes;
    d = [floeXArray(obsIdx(iObservation+1), obsFloesIdx).'; floeYArray(obsIdx(iObservation+1), obsFloesIdx).'; floeAlphaArray(obsIdx(iObservation+1), obsFloesIdx).'];
    E = mvnrnd(zeros(dObserved, 1), diag([sigobs^2*ones(1, nObsFloes), sigobs^2*ones(1, nObsFloes), sigalpha^2*ones(1, nObsFloes)]), nEnsemble)';
    
    D = repmat(d, 1, nEnsemble) + E;
    
    xObsArray = squeeze(xArray(obsFloesIdx, :));
    yObsArray = squeeze(yArray(obsFloesIdx, :));
    alphaObsArray = squeeze(alphaArray(obsFloesIdx, :));
    uObsArray = squeeze(uArray(obsFloesIdx, :));
    vObsArray = squeeze(vArray(obsFloesIdx, :));
    
    xMean = mean(xObsArray, 2);
    yMean = mean(yObsArray, 2);
    alphaMean = mean(alphaObsArray, 2);
    uMean = mean(xObsArray, 2);
    vMean = mean(yObsArray, 2);    
    
    Dp = D - [xObsArray; yObsArray; alphaObsArray];

    S = [xObsArray - xMean; yObsArray - yMean; alphaObsArray - alphaMean];

    C = S*S' + E*E';
    
    X = eye(nEnsemble) + S'/C * Dp;
    
    [Z, Lam] = eig(C);
    Cinv = Z * diag( 1./diag(Lam)) * Z.';
    X2 = diag( diag(Lam).^(-1/2) ) * Z.' * S;
    [U2, Sig2, V2] = svd(X2);
    
    oneN = ones(nEnsemble)/nEnsemble;
    
    lag = 10;
    delay = max(1, obsIdx(iObservation+1)-itsPerDay*lag);
    

    obsX = floeXArray(obsIdx(iObservation+1), obsFloesIdx).';
    obsY = floeYArray(obsIdx(iObservation+1), obsFloesIdx).';

    locRadFloe = 50*1e3;
    for iFloe = 1:nFloes
        if ~activeFloes(iFloe)
            continue
        end

         floeX = floeXArray(obsIdx(iObservation+1), iFloe);
         floeY = floeYArray(obsIdx(iObservation+1), iFloe);
         
         distX = obsX - floeX;
         distY = obsY - floeY;
         locIdx = distX.^2 + distY.^2 < locRadFloe^2;
         locObsIdx = locIdx;
         
         obsLocIdx = repmat(locObsIdx, nObsVar, 1);
         
         DpLoc = Dp(obsLocIdx, :);
         SLoc = S(obsLocIdx, :);
         ELoc = E(obsLocIdx, :);
         CLoc = SLoc*SLoc' + ELoc*ELoc';
         
         [ZLoc, LamLoc] = eig(CLoc);
         CinvLoc = ZLoc * diag( 1./diag(LamLoc)) * ZLoc.';
         X2Loc = diag( diag(LamLoc).^(-1/2) ) * ZLoc.' * SLoc;
         [U2Loc, Sig2Loc, V2Loc] = svd(X2Loc);
         oneN = ones(nEnsemble)/nEnsemble;
         XLoc = oneN + SLoc.' * CinvLoc * DpLoc * oneN + (eye(nEnsemble) - oneN) * V2Loc * sqrtm(eye(nEnsemble) - Sig2Loc.' * Sig2Loc) * V2Loc.';
  
         xEnsemble(iFloe, :, delay:end) = XLoc.' * squeeze(xEnsemble(iFloe, :, delay:end));
         yEnsemble(iFloe, :, delay:end) = XLoc.' * squeeze(yEnsemble(iFloe, :, delay:end));
         alphaEnsemble(iFloe, :, delay:end) = XLoc.' * squeeze(alphaEnsemble(iFloe, :, delay:end));
    
         xArray(iFloe, :) = XLoc.' * squeeze(xArray(iFloe, :)).';
         yArray(iFloe, :) = XLoc.' * squeeze(yArray(iFloe, :)).';
         
         alphaArray(iFloe, :) = XLoc.' * squeeze(alphaArray(iFloe, :)).';
         ksiArray(iFloe, :) = XLoc.' * squeeze(ksiArray(iFloe, :)).';
         uArray(iFloe, :) = XLoc.' * squeeze(uArray(iFloe, :)).';
         vArray(iFloe, :) = XLoc.' * squeeze(vArray(iFloe, :)).';
         dxArray(iFloe, :) = XLoc.' * squeeze(dxArray(iFloe, :)).';
         dyArray(iFloe, :) = XLoc.' * squeeze(dyArray(iFloe, :)).';
         dalphaArray(iFloe, :) = XLoc.' * squeeze(dalphaArray(iFloe, :)).';
         duArray(iFloe, :) = XLoc.' * squeeze(duArray(iFloe, :)).';
         dvArray(iFloe, :) = XLoc.' * squeeze(dvArray(iFloe, :)).';
         dksiArray(iFloe, :) = XLoc.' * squeeze(dksiArray(iFloe, :)).';
         
         thicknessArray(iFloe, :) = XLoc.' * squeeze(thicknessArray(iFloe, :)).';
         thicknessArray(iFloe, :) = max(thicknessArray(iFloe, :), 0.5);
         
        for iEnsemble = 1:nEnsemble
            floeEnsemble(iFloe, iEnsemble).Xi = xArray(iFloe, iEnsemble);
            floeEnsemble(iFloe, iEnsemble).Yi = yArray(iFloe, iEnsemble);
            
            floeEnsemble(iFloe, iEnsemble).alpha_i = alphaArray(iFloe, iEnsemble);
            floeEnsemble(iFloe, iEnsemble).Ui = uArray(iFloe, iEnsemble);
            floeEnsemble(iFloe, iEnsemble).Vi = vArray(iFloe, iEnsemble);
            floeEnsemble(iFloe, iEnsemble).ksi_ice = ksiArray(iFloe, iEnsemble);
            floeEnsemble(iFloe, iEnsemble).dXi_p = dxArray(iFloe, iEnsemble);
            floeEnsemble(iFloe, iEnsemble).dYi_p = dyArray(iFloe, iEnsemble);
            floeEnsemble(iFloe, iEnsemble).dalpha_i_p = dalphaArray(iFloe, iEnsemble);
            floeEnsemble(iFloe, iEnsemble).dUi_p = duArray(iFloe, iEnsemble);
            floeEnsemble(iFloe, iEnsemble).dVi_p = dvArray(iFloe, iEnsemble);
            floeEnsemble(iFloe, iEnsemble).dksi_ice_p = dksiArray(iFloe, iEnsemble);
            
            thickness = thicknessArray(iFloe, iEnsemble);
            floeEnsemble(iFloe, iEnsemble).mass = floeEnsemble(iFloe, iEnsemble).mass/floeEnsemble(iFloe, iEnsemble).thickness*thickness;
            floeEnsemble(iFloe, iEnsemble).inertia_moment = floeEnsemble(iFloe, iEnsemble).inertia_moment/floeEnsemble(iFloe, iEnsemble).thickness*thickness;
            floeEnsemble(iFloe, iEnsemble).thickness = thickness;
        end
    end
    
    thicknessEnsMean(:, iObservation+1) = mean(thicknessArray, 2);
    thicknessEnsStd(:, iObservation+1) = std(thicknessArray, 0, 2);


    % Update ocean
    pArray = zeros(nEnsModes, nEnsModes, nEnsemble);
    pInitialArray = zeros(nEnsModes, nEnsModes, nEnsemble);
    for iEnsemble = 1:nEnsemble
        pArray(:, :, iEnsemble) = ifft2(squeeze(phArray(:, :, iEnsemble)));
        pInitialArray(:, :, iEnsemble) = ifft2(squeeze(phInitial(:, :, iEnsemble)));
    end
    
    
    truncX = linspace(xMin*1e3, xMax*1e3, nEnsModes);
    truncY = linspace(yMin*1e3, yMax*1e3, nEnsModes);
   
    locRad = 200*1e3;
    for iMode = 1:nEnsModes
        for jMode = 1:nEnsModes
             gridX = truncX(jMode);
             gridY = truncY(iMode);
             
             distX = obsX - gridX;
             distY = obsY - gridY;
             locIdx = distX.^2 + distY.^2 < locRad^2;
             locObsIdx = locIdx;
             
             obsLocIdx = repmat(locObsIdx, nObsVar, 1);
             
             DpLoc = Dp(obsLocIdx, :);
             SLoc = S(obsLocIdx, :);
             ELoc = E(obsLocIdx, :);
             CLoc = SLoc*SLoc' + ELoc*ELoc';


             XLoc = eye(nEnsemble) + SLoc'/CLoc * DpLoc;
            pArray(iMode, jMode, :) = XLoc.' * squeeze(pArray(iMode, jMode, :));
            pInitialArray(iMode, jMode, :) = XLoc.' * squeeze(pInitialArray(iMode, jMode, :));
            

        end
    end
    for iEnsemble = 1:nEnsemble
        phArray(:, :, iEnsemble) = fft2(squeeze(pArray(:, :, iEnsemble)));
        phInitial(:, :, iEnsemble) = fft2(squeeze(pInitialArray(:, :, iEnsemble)));
    end

    phEnsMean(:, :, iObservation+1) = mean(phArray, 3);



    windArray = zeros(nWindModes, nWindModes, nEnsemble);
    windInitialArray = zeros(nWindModes, nWindModes, nEnsemble);
    for iEnsemble = 1:nEnsemble
        windArray(:, :, iEnsemble) = ifft2(squeeze(windHatArray(:, :, iEnsemble)));
        windInitialArray(:, :, iEnsemble) = ifft2(squeeze(windHatInitial(:, :, iEnsemble)));
    end
    

    truncX = linspace(xMin*1e3, xMax*1e3, nWindModes);
    truncY = linspace(yMin*1e3, yMax*1e3, nWindModes);

    for iMode = 1:nWindModes
        for jMode = 1:nWindModes

            XLoc = X;
            windArray(iMode, jMode, :) = XLoc.' * squeeze(windArray(iMode, jMode, :));
            windInitialArray(iMode, jMode, :) = XLoc.' * squeeze(windInitialArray(iMode, jMode, :));
            
        end
    end

    for iEnsemble = 1:nEnsemble
        windHatArray(:, :, iEnsemble) = fft2(squeeze(windArray(:, :, iEnsemble)));
        windHatInitial(:, :, iEnsemble) = fft2(squeeze(windInitialArray(:, :, iEnsemble)));
    end

    fprintf(' %d', iObservation)
end
fprintf('\n')







