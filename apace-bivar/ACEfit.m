function ACEfit_Par = ACEfit(ACEfit_Par)
%
% Element-wise parameter estimation by LR-SD and parameter inference
% with ERV estimator as the test statistic.
%

nMZ    = ACEfit_Par.nMZF*2;
nDZ    = ACEfit_Par.nDZF*2;
X      = ACEfit_Par.X;
Y      = ACEfit_Par.Y;
I_data = ACEfit_Par.I_data;

[nElm,n] = size(Y);

% Initalize result variables
[ERV,rhoG,rhoP] = deal(zeros(nElm));
[H2,Sig2]       = deal(zeros(nElm,1));
sdY             = zeros(nElm,nElm,4);

[ZTZ,pinvZTZ]   = Mk_XTXpinv(nMZ/2,nDZ/2,0,n);

pinvXTX = X*pinv(X);
p       = size(X,2);


%
% Heritability estimation
%
for i = I_data
    
    iY    = Y(i,:)';
    iRes  = iY - pinvXTX*iY; 
    iSig2 = sum(iRes.^2)/(n-p);
    % iSig2 = var(iRes);
    
    [hh,~]  = LRSD_ACE_FAM([1:2:nMZ; 2:2:nMZ]',nMZ+[1:2:nDZ; 2:2:nDZ]',zeros(0,2),iRes,ZTZ,pinvZTZ,iSig2);
    
    Sig2(i) = iSig2;
    H2(i)   = hh(1)/sum(hh);
    
end

%
% Parameter estimation
%
for id = 1:length(I_data)
    
    i    = I_data(id); 
    
    iY   = Y(i,:)';
    iRes = iY - pinvXTX*iY;
    iStd = iRes/std(iRes);
    % iStd = iRes/sqrt(Sig2(i));
    
    for j = I_data(id+1:end)
        
        jY   = Y(j,:)';
        jRes = jY - pinvXTX*jY;
        jStd = jRes/std(jRes);
                
        [sd,rhoP12,ERV12] = LRSD_BIVAR(nMZ,nDZ,n,iStd,jStd,1,1);
        % [sd,rhoP12,ERV12] = LRSD_BIVAR(nMZ,nDZ,n,iStd,jStd,var(iStd),var(jStd));
        
        sdY(i,j,:) = sd;
        
        ERV(i,j)  = abs(ERV12);
        rhoP(i,j) = rhoP12;
        if ( H2(i)>0 && H2(j)>0 )
            rhoG(i,j) = min(max(ERV12/sqrt(H2(i)*H2(j)),-1),1);
        end
        
    end
    
end

% Save for later computation
ACEfit_Par.sdY    = sdY;

% Obtain the test statistic: ERV
ACEfit_Par.Stats = ERV;
ACEfit_Par.mT    = max(max(ACEfit_Par.Stats));

% Save result file
save(fullfile(ACEfit_Par.ResDir,'Results'),'H2','Sig2','rhoP','rhoG','ERV','-v7.3');

fprintf('.\n');

return
