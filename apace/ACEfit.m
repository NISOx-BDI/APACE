function ACEfit_Par = ACEfit(ACEfit_Par)
%
% Element- or voxel-wise phenotypic heritability estimation using LR-SD and
% heritability inference with likelihood ratio test (LRT).
%
% The LRT statistic is chosen as the test statistic rather than h2
% estimator as theoretically h2 estimator is not sufficient nor pivotal
% while LRT statistic is pivotal. The distribution of LRT statistic does
% not depend on unknown nuisance parameters and remains the same under the
% null hypothesis, and thus the pivotal test statistic is recommendated in
% permutation test.
%

X = ACEfit_Par.X;
Y = ACEfit_Par.Y;

nFam = ACEfit_Par.nFam;
nMZF = ACEfit_Par.nMZF;

I_MZ   = ACEfit_Par.I_MZ;
I_DZ   = ACEfit_Par.I_DZ;
I_Sib  = ACEfit_Par.I_Sib;
I_data = ACEfit_Par.I_data;

[nElm,n] = size(Y);

% Make residuals, compute total variance, initialise
Res           = Y - Y*(X*pinv(X))';
Sigma2        = sum(Res.^2,2)/(n-size(X,2));
[XTX,pinvXTX] = Mk_XTXpinv(size(I_MZ,1),size(I_DZ,1),size(I_Sib,1),n);

BadElm        = 0;

switch upper(ACEfit_Par.Model)
    
    case 'ACE'
        
        % Initalize results variables
        [h2,c2,e2] = deal(zeros(nElm,1));
        if ~ACEfit_Par.NoImg
            [Stats,Stats_C]       = deal(zeros(nElm,1));
            [AsyPval_A,AsyPval_C] = deal(ones(nElm,1));
        end
        
        for i = I_data
            
            iY   = Y(i,:)';
            iRes = Res(i,:)';
            sig2 = Sigma2(i);
            
            %
            % Heritability estimation
            %
            [hha,HH] = LRSD_ACE_FAM(I_MZ,I_DZ,I_Sib,iRes,XTX,pinvXTX,sig2);
            hha      = hha/sum(hha);
            
            h2(i)    = hha(1);
            c2(i)    = hha(2);
            e2(i)    = hha(3);
            
            if ~ACEfit_Par.NoImg
                
                %
                % Likelihood ratio test
                %
                
                if (hha(3)<1e-5)
                    % Zero environmental variance... give up
                    BadElm       = BadElm + 1;
                    AsyPval_A(i) = 1;
                    AsyPval_C(i) = 1;
                    Stats(i)     = 0;
                    Stats_C(i)   = 0;
                else
                    
                    if ( min(HH(:,3))>=0 )
                        hh0 = HH(:,3);
                    else
                        hh0 = HH(:,1);
                    end
                    if ( min(HH(:,2))>=0 )
                        hh_AE = HH(:,2);
                    else
                        hh_AE = HH(:,1);
                    end
                    hh0   = hh0/sum(hh0);
                    hh_AE = hh_AE/sum(hh_AE);
                    
                    [det_Va,det_V0,det_V_AE] = deal(1);
                    [invMZR,invOTR,invR0,invMZR_AE,invOTR_AE] = deal(cell(max(nFam),1));
                    [invMZR{1},invOTR{1},invR0{1},invMZR_AE{1},invOTR_AE{1}] = deal(1);
                    
                    for nF=2:max(nFam)
                        n_nF          = sum(nFam==nF);
                        n_nFMZ        = sum(nFam(1:nMZF)==nF);
                        
                        %%% ACE model
                        a             = hha(1)/2+hha(3);
                        b             = hha(1)/2+hha(2);
                        R             = eye(nF)*a + ones(nF)*b;
                        det_Va        = det_Va*det(R)^(n_nF-n_nFMZ);
                        invOTR{nF}    = (eye(nF) - ones(nF)*b/(a+b*nF))/a;
                        
                        R(1,2)        = hha(1)+hha(2);
                        R(2,1)        = R(1,2);
                        det_Va        = det_Va*det(R)^n_nFMZ;
                        invMZR{nF}    = sparse(pinv(R));
                        
                        %%% CE model
                        R0            = eye(nF)*hh0(3) + ones(nF)*hh0(2);
                        det_V0        = det_V0*det(R0)^n_nF;
                        invR0{nF}     = (eye(nF) - ones(nF)*hh0(2)/(hh0(3)+hh0(2)*nF))/hh0(3);
                        
                        %%% AE model
                        a             = hh_AE(1)/2+hh_AE(3);
                        b             = hh_AE(1)/2;
                        R_AE          = eye(nF)*a + ones(nF)*b;
                        det_V_AE      = det_V_AE*det(R_AE)^(n_nF-n_nFMZ);
                        invOTR_AE{nF} = (eye(nF) - ones(nF)*b/(a+b*nF))/a;
                        
                        R_AE(1,2)     = hh_AE(1);
                        R_AE(2,1)     = R_AE(1,2);
                        det_V_AE      = det_V_AE*det(R_AE)^n_nFMZ;
                        invMZR_AE{nF} = sparse(pinv(R_AE));
                    end
                    
                    % ACE model
                    Ua  = blkdiag(invMZR{nFam(1:nMZF)},invOTR{nFam(nMZF+1:end)});
                    Za  = pinv(X'*Ua*X);
                    Pa  = Ua-Ua*X*Za*X'*Ua;
                    rla = -(iY'*Pa*iY/sig2+log(det_Va)+log(det(X'*Ua*X)))/2;
                    
                    % CE model
                    U0  = blkdiag(invR0{nFam});
                    Z0  = pinv(X'*U0*X);
                    P0  = U0-U0*X*Z0*X'*U0;
                    rl0 = -(iY'*P0*iY/sig2+log(det_V0)+log(det(X'*U0*X)))/2;
                    
                    % AE model
                    U_AE  = blkdiag(invMZR_AE{nFam(1:nMZF)}, invOTR_AE{nFam(nMZF+1:end)});
                    Z_AE  = pinv(X'*U_AE*X);
                    P_AE  = U_AE-U_AE*X*Z_AE*X'*U_AE;
                    rl_AE = -(iY'*P_AE*iY/sig2+log(det_V_AE)+log(det(X'*U_AE*X)))/2;
                    
                    % Compute LRT test statistic
                    T   = -2*(rl0-rla);
                    T_C = -2*(rl_AE-rla);
                    
                    % Obtain the resulting asymptotic p-value
                    if T>0
                        AsyPval_A(i) = (1-chi2cdf(T,1))/2;
                    end
                    if T_C>0
                        AsyPval_C(i) = (1-chi2cdf(T_C,1))/2;
                    end
                    
                    % Save LRT test statistic
                    Stats(i)   = max(0,T);
                    Stats_C(i) = max(0,T_C);
                    
                end
            end
            
        end
        
        
        h2O      = h2(I_data);
        c2O      = c2(I_data);
        e2O      = e2(I_data);
        Sigma2O  = Sigma2(I_data);
        
        AO       = Sigma2O.*h2O;
        wh2O     = mean(AO)/mean(Sigma2O);
        
        quartile = quantile(h2O,[0.50 0.75]);
        medh2O   = quartile(1);
        q3h2O    = quartile(2);
        meanh2O  = mean(h2O);
        mGmedh2O = mean(h2O(h2O>medh2O));
        % if isnan(mGmedh2O)
        %     mGmedh2O = mean(h2O(h2O>=medh2O));
        % end
        mGq3h2O  = mean(h2O(h2O>q3h2O));
        % if isnan(mGq3h2O)
        %     mGq3h2O  = mean(h2O(h2O>=q3h2O));
        % end
        
        % Summary statistics: meanh2, wh2, median, q3, mean(h2>median), mean(h2>q3)
        ACEfit_Par.SummaryA = [meanh2O; wh2O; medh2O; q3h2O; mGmedh2O; mGq3h2O];
        
        CO       = Sigma2O.*c2O;
        wc2O     = mean(CO)/mean(Sigma2O);
        
        quartile = quantile(c2O,[0.50 0.75]);
        medc2O   = quartile(1);
        q3c2O    = quartile(2);
        meanc2O  = mean(c2O);
        mGmedc2O = mean(c2O(c2O>medc2O));
        % if isnan(mGmedc2O)
        %     mGmedc2O = mean(c2O(c2O>=medc2O));
        % end
        mGq3c2O  = mean(c2O(c2O>q3c2O));
        % if isnan(mGq3c2O)
        %     mGq3c2O  = mean(c2O(c2O>=q3c2O));
        % end
        
        % Summary statistics: meanc2, wc2, median, q3, mean(c2>median), mean(c2>q3)
        ACEfit_Par.SummaryC = [meanc2O; wc2O; medc2O; q3c2O; mGmedc2O; mGq3c2O];
        
        EO       = Sigma2O.*e2O;
        we2O     = mean(EO)/mean(Sigma2O);
        
        quartile = quantile(e2O,[0.50 0.75]);
        mede2O   = quartile(1);
        q3e2O    = quartile(2);
        meane2O  = mean(e2O);
        mGmede2O = mean(e2O(e2O>mede2O));
        % if isnan(mGmede2O)
        %     mGmede2O = mean(e2O(e2O>=mede2O));
        % end
        mGq3e2O  = mean(e2O(e2O>q3e2O));
        % if isnan(mGq3e2O)
        %     mGq3e2O  = mean(e2O(e2O>=q3e2O));
        % end
        
        
        % Summary statistics: meane2, we2, median, q3, mean(e2>median), mean(e2>q3)
        ACEfit_Par.SummaryE = [meane2O; we2O; mede2O; q3e2O; mGmede2O; mGq3e2O];
        
        % Write out the ouptut images
        WriteData(h2,           ACEfit_Par.Vs, 'ACE_A_h2', ACEfit_Par.ResDir);
        WriteData(c2,           ACEfit_Par.Vs, 'ACE_C_c2', ACEfit_Par.ResDir);
        WriteData(e2,           ACEfit_Par.Vs, 'ACE_E_e2', ACEfit_Par.ResDir);
        WriteData(sqrt(Sigma2), ACEfit_Par.Vs, 'Stdev',    ACEfit_Par.ResDir);
        
        
    case 'AE'
        
        % Initalize results variables
        [h2,e2] = deal(zeros(nElm,1));
        if ~ACEfit_Par.NoImg
            Stats     = zeros(nElm,1);
            AsyPval_A = ones(nElm,1);
        end
        
        P0  = eye(n)-X*pinv(X);
        
        for i = I_data
            
            iY   = Y(i,:)';
            iRes = Res(i,:)';
            sig2 = Sigma2(i);
            
            %
            % Heritability estimation
            %
            [~,HH] = LRSD_ACE_FAM(I_MZ,I_DZ,I_Sib,iRes,XTX,pinvXTX,sig2);
            if ( min(HH(:,2))>=0 )
                hha = HH(:,2);
            else
                hha = HH(:,1);
            end
            hha = hha/sum(hha);
            
            h2(i) = hha(1);
            e2(i) = hha(3);
            
            if ~ACEfit_Par.NoImg
                
                %
                % Likelihood ratio test
                %
                
                if (hha(3)<1e-5)
                    % Zero enviornmental variance... give up
                    BadElm       = BadElm + 1;
                    AsyPval_A(i) = 1;
                    Stats(i)     = 0;
                else
                    
                    index = ACEcode(hha);
                    if ( index==ACEcode([0 0 1]) )
                        % Parameter A is zero, test-statistic is zero
                        T = 0;
                    else
                        det_Va                = 1;
                        [invMZR,invOTR]       = deal(cell(max(nFam),1));
                        [invMZR{1},invOTR{1}] = deal(1);
                        
                        for nF=2:max(nFam)
                            n_nF       = sum(nFam==nF);
                            n_nFMZ     = sum(nFam(1:nMZF)==nF);
                            
                            %%% AE model
                            a          = hha(1)/2+hha(3);
                            b          = hha(1)/2;
                            Ra         = eye(nF)*a + ones(nF)*b;
                            det_Va     = det_Va*det(Ra)^(n_nF-n_nFMZ);
                            invOTR{nF} = (eye(nF) - ones(nF)*b/(a+b*nF))/a;
                            
                            Ra(1,2)    = hha(1);
                            Ra(2,1)    = Ra(1,2);
                            det_Va     = det_Va*det(Ra)^n_nFMZ;
                            invMZR{nF} = sparse(pinv(Ra));
                        end
                        
                        % AE model
                        Ua  = blkdiag(invMZR{nFam(1:nMZF)}, invOTR{nFam(nMZF+1:end)});
                        Za  = pinv(X'*Ua*X);
                        Pa  = Ua-Ua*X*Za*X'*Ua;
                        rla = -(iY'*Pa*iY/sig2+log(det_Va)+log(det(X'*Ua*X)))/2;
                        
                        % E model
                        rl0 = -(iY'*P0*iY/sig2+log(det(X'*X)))/2;
                        
                        % Compute LRT test statistic
                        T   = -2*(rl0-rla);
                    end
                    
                    % Obtain the resulting asymptotic p-value
                    if T>0
                        AsyPval_A(i) = (1-chi2cdf(T,1))/2;
                    end
                    
                    % Save LRT test statistic
                    Stats(i) = max(0,T);
                    
                end
            end
            
        end
        
        
        h2O      = h2(I_data);
        e2O      = e2(I_data);
        Sigma2O  = Sigma2(I_data);
        
        AO       = Sigma2O.*h2O;
        wh2O     = mean(AO)/mean(Sigma2O);
        
        quartile = quantile(h2O,[0.50 0.75]);
        medh2O   = quartile(1);
        q3h2O    = quartile(2);
        meanh2O  = mean(h2O);
        mGmedh2O = mean(h2O(h2O>medh2O));
        % if isnan(mGmedh2O)
        %     mGmedh2O = mean(h2O(h2O>=medh2O));
        % end
        mGq3h2O  = mean(h2O(h2O>q3h2O));
        % if isnan(mGq3h2O)
        %     mGq3h2O  = mean(h2O(h2O>=q3h2O));
        % end
        
        % Summary statistics: meanh2, wh2, median, q3, mean(h2>median), mean(h2>q3)
        ACEfit_Par.SummaryA = [meanh2O; wh2O; medh2O; q3h2O; mGmedh2O; mGq3h2O];
        
        EO       = Sigma2O.*e2O;
        we2O     = mean(EO)/mean(Sigma2O);
        
        quartile = quantile(e2O,[0.50 0.75]);
        mede2O   = quartile(1);
        q3e2O    = quartile(2);
        meane2O  = mean(e2O);
        mGmede2O = mean(e2O(e2O>mede2O));
        % if isnan(mGmede2O)
        %     mGmede2O = mean(e2O(e2O>=mede2O));
        % end
        mGq3e2O  = mean(e2O(e2O>q3e2O));
        % if isnan(mGq3e2O)
        %     mGq3e2O  = mean(e2O(e2O>=q3e2O));
        % end
        
        % Summary statistics: meane2, we2, median, q3, mean(e2>median), mean(e2>q3)
        ACEfit_Par.SummaryE = [meane2O; we2O; mede2O; q3e2O; mGmede2O; mGq3e2O];
        
        % Write out the ouptut images
        WriteData(h2,           ACEfit_Par.Vs, 'ACE_A_h2', ACEfit_Par.ResDir);
        WriteData(e2,           ACEfit_Par.Vs, 'ACE_E_e2', ACEfit_Par.ResDir);
        WriteData(sqrt(Sigma2), ACEfit_Par.Vs, 'Stdev',    ACEfit_Par.ResDir);
end

if BadElm>0
    fprintf('WARNING: Voxels/elements found with zero environmental variance (%d).  Check for serverely sparse discrete data.\n',BadElm);
end

% Plot the distribution of overall h2 estimates
f = 0;
% 1) h2_dist
f = f+1;
figure(f);
plot(sort(h2O));
grid on
title('Heritabilty h^2')
axis tight
set(gcf,'PaperPosition',[0 0 10 6])
set(gcf,'PaperSize',[10 6])
print('-dpdf',fullfile(ACEfit_Par.ResDir,'h2_dist.pdf'));
% 2) h2_hist
f = f+1;
figure(f);
hist(h2O(h2O>0),50);
per = mean(h2O>0);
title('Heritabilty h^2 (h^2 > 0 only)')
xlabel(sprintf('Positive h^2 (%d%% of h^2 > 0)',round(per*100)))
axis tight
set(gcf,'PaperPosition',[0 0 10 6])
set(gcf,'PaperSize',[10 6])
print('-dpdf',fullfile(ACEfit_Par.ResDir,'h2_hist.pdf'));



%
% Image-wise inference
%
if ~ACEfit_Par.NoImg
    
    switch upper(ACEfit_Par.Model)
        case 'ACE'
            % Save maximum LRT statistic (mT)
            ACEfit_Par.mT    = max(Stats);
            ACEfit_Par.Stats = Stats;
            WriteData(Stats,             ACEfit_Par.Vs, 'ACE_A_LRT',           ACEfit_Par.ResDir);
            WriteData(Stats_C,           ACEfit_Par.Vs, 'ACE_C_LRT',           ACEfit_Par.ResDir);
            WriteData(-log10(AsyPval_A), ACEfit_Par.Vs, 'ACE_A_LRT_vox_Pasym', ACEfit_Par.ResDir);
            WriteData(-log10(AsyPval_C), ACEfit_Par.Vs, 'ACE_C_LRT_vox_Pasym', ACEfit_Par.ResDir);
        case 'AE'
            % Save maximum LRT statistic (mT)
            ACEfit_Par.mT    = max(Stats);
            ACEfit_Par.Stats = Stats;
            WriteData(Stats,             ACEfit_Par.Vs, 'ACE_A_LRT',           ACEfit_Par.ResDir);
            WriteData(-log10(AsyPval_A), ACEfit_Par.Vs, 'ACE_A_LRT_vox_Pasym', ACEfit_Par.ResDir);
    end
    
    
    %
    % Cluster inference
    %
    if ACEfit_Par.Vs.ClustInf
        
        if ( ~isfield(ACEfit_Par,'alpha_CFT') || isempty(ACEfit_Par.alpha_CFT) )
            ACEfit_Par.alpha_CFT = 0.05;
        end
        
        if length(ACEfit_Par.Vs.Dim)==1
            Dim = [ACEfit_Par.Vs.Dim 1];
        else
            Dim = ACEfit_Par.Vs.Dim;
        end
        
        Tclus  = spm_invXcdf(1-2*ACEfit_Par.alpha_CFT,1);
        Tstats = reshape(Stats,Dim);
        Tstats(Tstats<Tclus) = 0;
        
        % Cluster-forming procedure
        [L,NUM] = spm_bwlabel(Tstats,18);
        
        if NUM>0
            cluster_size = zeros(NUM,1);
            cluster_mass = zeros(NUM,1);
            for i=1:NUM
                cluster_size(i) = sum(L(:)==i);
                cluster_mass(i) = sum(Tstats(L(:)==i));
            end
            mK = max(cluster_size);
            mM = max(cluster_mass);
        else
            mK = 0;
            mM = 0;
        end
        % Save maximum cluster statistics of size (mK) and mass (mM)
        ACEfit_Par.mK = mK;
        ACEfit_Par.mM = mM;
        
    end
    
end

return
