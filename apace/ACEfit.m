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
%_______________________________________________________________________
% Version: http://github.com/NISOx-BDI/APACE/tree/$Format:%h$
%          $Format:%ci$

X = ACEfit_Par.X;
Y = ACEfit_Par.Y;

nFam = ACEfit_Par.nFam;
nMZF = ACEfit_Par.nMZF;

I_MZ   = ACEfit_Par.I_MZ;
I_DZ   = ACEfit_Par.I_DZ;
I_Sib  = ACEfit_Par.I_Sib;
I_data = ACEfit_Par.I_data;

[nElm,n] = size(Y);

% Initalize results variables
[hh,hhAE,hhCE] = deal(zeros(nElm,3));

% Make residuals, compute total variance, initialise
Res    = Y - Y*(X*pinv(X))';
Sigma2 = sum(Res.^2,2)/(n-size(X,2));

BadElm = 0;

%
% Heritability estimation
%
[hhACE0,hhAE0,hhCE0] = LRSD_ACE_VEC(I_MZ,I_DZ,I_Sib,Res(I_data,:),Sigma2(I_data));

switch upper(ACEfit_Par.Model)

    case 'ACE'

        hh(I_data,:) = hhACE0./repmat(sum(hhACE0,2),1,3);

        h2 = hh(:,1);
        c2 = hh(:,2);
        e2 = hh(:,3);

        %
        % Image-wise inference
        %
        if ~ACEfit_Par.NoImg

            %
            % Likelihood ratio test
            %
            hhCE(I_data,:) = hhCE0./repmat(sum(hhCE0,2),1,3);
            hhAE(I_data,:) = hhAE0./repmat(sum(hhAE0,2),1,3);

            % Initialization
            [T,T_C]               = deal(zeros(nElm,1));
            [AsyPval_A,AsyPval_C] = deal(ones(nElm,1));

            for i = I_data

                iY    = Y(i,:)';
                hha   = hh(i,:)';
                hh0   = hhCE(i,:)';
                hh_AE = hhAE(i,:)';
                sig2  = Sigma2(i);

                if (hha(3)<1e-5)
                    % Zero environmental variance... give up
                    BadElm       = BadElm + 1;
                    AsyPval_A(i) = 1;
                    AsyPval_C(i) = 1;
                    T(i)         = 0;
                    T_C(i)       = 0;
                else
                    index = ACEcode(hha);
                    if ( index~=ACEcode([0 0 1]) )  % Found model is not E

                        [det_Va,det_V0,det_V_AE] = deal(1);
                        [invMZR,invOTR,invR0,invMZR_AE,invOTR_AE] = deal(cell(max(nFam),1));
                        [invMZR{1},invOTR{1},invR0{1},invMZR_AE{1},invOTR_AE{1}] = deal(1);

                        for nF = 2:max(nFam)
                            n_nF                  = sum(nFam==nF);
                            n_nFMZ                = sum(nFam(1:nMZF)==nF);

                            %%% ACE model
                            a                     = hha(1)/2+hha(3);
                            b                     = hha(1)/2+hha(2);
                            R                     = eye(nF)*a + ones(nF)*b;
                            det_Va                = det_Va*det(R)^(n_nF-n_nFMZ);
                            invOTR{nF}            = (eye(nF) - ones(nF)*b/(a+b*nF))/a;

                            [R(1,2),R(2,1)]       = deal(hha(1)+hha(2));
                            det_Va                = det_Va*det(R)^n_nFMZ;
                            invMZR{nF}            = sparse(pinv(R));

                            %%% CE model
                            R0                    = eye(nF)*hh0(3) + ones(nF)*hh0(2);
                            det_V0                = det_V0*det(R0)^n_nF;
                            invR0{nF}             = (eye(nF) - ones(nF)*hh0(2)/(hh0(3)+hh0(2)*nF))/hh0(3);

                            %%% AE model
                            a                     = hh_AE(1)/2+hh_AE(3);
                            b                     = hh_AE(1)/2;
                            R_AE                  = eye(nF)*a + ones(nF)*b;
                            det_V_AE              = det_V_AE*det(R_AE)^(n_nF-n_nFMZ);
                            invOTR_AE{nF}         = (eye(nF) - ones(nF)*b/(a+b*nF))/a;

                            [R_AE(1,2),R_AE(2,1)] = deal(hh_AE(1));
                            det_V_AE              = det_V_AE*det(R_AE)^n_nFMZ;
                            invMZR_AE{nF}         = sparse(pinv(R_AE));
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
                        T(i)   = -2*(rl0-rla);
                        T_C(i) = -2*(rl_AE-rla);
                    end

                    if ( index==ACEcode([0 1 1]) ) % Found CE model
                        T(i)   = 0;
                    elseif ( index==ACEcode([1 0 1]) ) % Found AE model
                        T_C(i) = 0;
                    end

                end
            end

            % Save LRT test statistic
            Stats   = max(0,T);
            Stats_C = max(0,T_C);

            % Obtain the resulting asymptotic p-value
            AsyPval_A(T>0) = (1-chi2cdf(T(T>0),1))/2;
            AsyPval_C(T_C>0) = (1-chi2cdf(T_C(T_C>0),1))/2;

            % Write out the LRT statistic images
            WriteData(Stats,             ACEfit_Par.Vs, 'ACE_A_LRT',           ACEfit_Par.ResDir);
            WriteData(Stats_C,           ACEfit_Par.Vs, 'ACE_C_LRT',           ACEfit_Par.ResDir);
            WriteData(-log10(AsyPval_A), ACEfit_Par.Vs, 'ACE_A_LRT_vox_Pasym', ACEfit_Par.ResDir);
            WriteData(-log10(AsyPval_C), ACEfit_Par.Vs, 'ACE_C_LRT_vox_Pasym', ACEfit_Par.ResDir);

            % Save maximum LRT statistic (mT)
            ACEfit_Par.mT    = max(Stats);
            ACEfit_Par.Stats = Stats;

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

        %%% Compute the summary statistics
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
        mGmedh2O = mean(h2O(h2O>=medh2O));
        mGq3h2O  = mean(h2O(h2O>=q3h2O));

        % Summary statistics: meanh2, wh2, median, q3, mean(h2>median), mean(h2>q3)
        ACEfit_Par.SummaryA = [meanh2O; wh2O; medh2O; q3h2O; mGmedh2O; mGq3h2O];

        CO       = Sigma2O.*c2O;
        wc2O     = mean(CO)/mean(Sigma2O);

        quartile = quantile(c2O,[0.50 0.75]);
        medc2O   = quartile(1);
        q3c2O    = quartile(2);
        meanc2O  = mean(c2O);
        mGmedc2O = mean(c2O(c2O>=medc2O));
        mGq3c2O  = mean(c2O(c2O>=q3c2O));

        % Summary statistics: meanc2, wc2, median, q3, mean(c2>median), mean(c2>q3)
        ACEfit_Par.SummaryC = [meanc2O; wc2O; medc2O; q3c2O; mGmedc2O; mGq3c2O];

        EO       = Sigma2O.*e2O;
        we2O     = mean(EO)/mean(Sigma2O);

        quartile = quantile(e2O,[0.50 0.75]);
        mede2O   = quartile(1);
        q3e2O    = quartile(2);
        meane2O  = mean(e2O);
        mGmede2O = mean(e2O(e2O>=mede2O));
        mGq3e2O  = mean(e2O(e2O>=q3e2O));

        % Summary statistics: meane2, we2, median, q3, mean(e2>median), mean(e2>q3)
        ACEfit_Par.SummaryE = [meane2O; we2O; mede2O; q3e2O; mGmede2O; mGq3e2O];

        % Write out the variance parameter estimate images
        WriteData(h2,           ACEfit_Par.Vs, 'ACE_A_h2', ACEfit_Par.ResDir);
        WriteData(c2,           ACEfit_Par.Vs, 'ACE_C_c2', ACEfit_Par.ResDir);
        WriteData(e2,           ACEfit_Par.Vs, 'ACE_E_e2', ACEfit_Par.ResDir);
        WriteData(sqrt(Sigma2), ACEfit_Par.Vs, 'Stdev',    ACEfit_Par.ResDir);

    case 'AE'

        hh(I_data,:) = hhAE0./repmat(sum(hhAE0,2),1,3);

        h2 = hh(:,1);
        e2 = hh(:,3);

        %
        % Image-wise inference
        %
        if ~ACEfit_Par.NoImg

            %
            % Likelihood ratio test
            %

            % Initialization
            T         = zeros(nElm,1);
            AsyPval_A = ones(nElm,1);

            P0        = eye(n)-X*pinv(X);

            for i = I_data

                iY   = Y(i,:)';
                hha  = hh(i,:)';
                sig2 = Sigma2(i);

                if (hha(3)<1e-5)
                    % Zero enviornmental variance... give up
                    BadElm       = BadElm + 1;
                    AsyPval_A(i) = 1;
                    T(i)         = 0;
                else

                    % Parameter A is non-zero... compute log-likelihood...
                    index = ACEcode(hha);
                    if ( index==ACEcode([1 0 1]) )

                        det_Va                = 1;
                        [invMZR,invOTR]       = deal(cell(max(nFam),1));
                        [invMZR{1},invOTR{1}] = deal(1);

                        for nF = 2:max(nFam)
                            n_nF              = sum(nFam==nF);
                            n_nFMZ            = sum(nFam(1:nMZF)==nF);

                            %%% AE model
                            a                 = hha(1)/2+hha(3);
                            b                 = hha(1)/2;
                            Ra                = eye(nF)*a + ones(nF)*b;
                            det_Va            = det_Va*det(Ra)^(n_nF-n_nFMZ);
                            invOTR{nF}        = (eye(nF) - ones(nF)*b/(a+b*nF))/a;

                            [Ra(1,2),Ra(2,1)] = deal(hha(1));
                            det_Va            = det_Va*det(Ra)^n_nFMZ;
                            invMZR{nF}        = sparse(pinv(Ra));
                        end

                        % AE model
                        Ua  = blkdiag(invMZR{nFam(1:nMZF)}, invOTR{nFam(nMZF+1:end)});
                        Za  = pinv(X'*Ua*X);
                        Pa  = Ua-Ua*X*Za*X'*Ua;
                        rla = -(iY'*Pa*iY/sig2+log(det_Va)+log(det(X'*Ua*X)))/2;

                        % E model
                        rl0 = -(iY'*P0*iY/sig2+log(det(X'*X)))/2;

                        % Compute LRT test statistic
                        T(i) = -2*(rl0-rla);
                    end

                end
            end

            % Save LRT test statistic
            Stats = max(0,T);

            % Obtain the resulting asymptotic p-value
            AsyPval_A(T>0) = (1-chi2cdf(T(T>0),1))/2;

            % Write out the LRT statistic images
            WriteData(Stats,             ACEfit_Par.Vs, 'AE_A_LRT',           ACEfit_Par.ResDir);
            WriteData(-log10(AsyPval_A), ACEfit_Par.Vs, 'AE_A_LRT_vox_Pasym', ACEfit_Par.ResDir);

            % Save maximum LRT statistic (mT)
            ACEfit_Par.mT    = max(Stats);
            ACEfit_Par.Stats = Stats;

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

        %%% Compute the summary statistics
        h2O      = h2(I_data);
        e2O      = e2(I_data);
        Sigma2O  = Sigma2(I_data);

        AO       = Sigma2O.*h2O;
        wh2O     = mean(AO)/mean(Sigma2O);

        quartile = quantile(h2O,[0.50 0.75]);
        medh2O   = quartile(1);
        q3h2O    = quartile(2);
        meanh2O  = mean(h2O);
        mGmedh2O = mean(h2O(h2O>=medh2O));
        mGq3h2O  = mean(h2O(h2O>=q3h2O));

        % Summary statistics: meanh2, wh2, median, q3, mean(h2>median), mean(h2>q3)
        ACEfit_Par.SummaryA = [meanh2O; wh2O; medh2O; q3h2O; mGmedh2O; mGq3h2O];

        EO       = Sigma2O.*e2O;
        we2O     = mean(EO)/mean(Sigma2O);

        quartile = quantile(e2O,[0.50 0.75]);
        mede2O   = quartile(1);
        q3e2O    = quartile(2);
        meane2O  = mean(e2O);
        mGmede2O = mean(e2O(e2O>=mede2O));
        mGq3e2O  = mean(e2O(e2O>=q3e2O));

        % Summary statistics: meane2, we2, median, q3, mean(e2>median), mean(e2>q3)
        ACEfit_Par.SummaryE = [meane2O; we2O; mede2O; q3e2O; mGmede2O; mGq3e2O];

        % Write out the variance parameter estimate images
        WriteData(h2,           ACEfit_Par.Vs, 'AE_A_h2', ACEfit_Par.ResDir);
        WriteData(e2,           ACEfit_Par.Vs, 'AE_E_e2', ACEfit_Par.ResDir);
        WriteData(sqrt(Sigma2), ACEfit_Par.Vs, 'Stdev',   ACEfit_Par.ResDir);

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
title('Heritabilty h^2 - Cumulative Distribution')
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


return
