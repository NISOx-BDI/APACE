function [mK,mM,mT,uCount,SummaryA] = ACEfit_Perm(ACEfit_Par,Perm_label)
%
% Permutations
%
% Perm_label - Permutation relabeling
%

X  = ACEfit_Par.X;
Y  = ACEfit_Par.Y;
Vs = ACEfit_Par.Vs;

I_MZ   = ACEfit_Par.I_MZ;
I_DZ   = ACEfit_Par.I_DZ;
I_Sib  = ACEfit_Par.I_Sib;
I_data = ACEfit_Par.I_data;

nFam = ACEfit_Par.nFam;
nMZF = ACEfit_Par.nMZF;

[nElm,n] = size(Y);

% Permutation relabeling
I_TWIN = [I_MZ; I_DZ];
I_TWIN = I_TWIN(Perm_label,:);
I_MZ   = I_TWIN(1:nMZF,:);
I_DZ   = I_TWIN((nMZF+1):end,:);

% Initalize results variables
[mK,mM,mT,uCount] = deal(0);

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
        
        [hh,hhCE]    = deal(zeros(nElm,3));
        hh(I_data,:) = hhACE0./repmat(sum(hhACE0,2),1,3);
        
        %
        % Image-wise inference
        %
        if ~ACEfit_Par.NoImg
            
            %
            % Likelihood ratio test
            %
            hhCE(I_data,:) = hhCE0./repmat(sum(hhCE0,2),1,3);
            
            % Initialization
            Stats = zeros(nElm,1);
            
            for i = I_data
                
                iY   = Y(i,:)';
                hha  = hh(i,:)';
                hh0  = hhCE(i,:)';
                sig2 = Sigma2(i);
                
                if (hha(3)<1e-5)
                    % Zero environmental variance... give up
                    BadElm   = BadElm + 1;
                    Stats(i) = 0;
                else
                    
                    % Parameter A is non-zero... compute log-likelihood...
                    index = ACEcode(hha);
                    if ( index==ACEcode([1 0 1]) || index==ACEcode([1 1 1]) )
                        
                        [det_Va,det_V0] = deal(1);
                        [invMZR,invOTR,invR0] = deal(cell(max(nFam),1));
                        [invMZR{1},invOTR{1},invR0{1}] = deal(1);
                        
                        for nF = 2:max(nFam)
                            n_nF   = sum(nFam==nF);
                            n_nFMZ = sum(nFam(Perm_label(1:nMZF))==nF);
                            
                            %%% ACE model
                            a               = hha(1)/2+hha(3);
                            b               = hha(1)/2+hha(2);
                            R               = eye(nF)*a + ones(nF)*b;
                            det_Va          = det_Va*det(R)^(n_nF-n_nFMZ);
                            invOTR{nF}      = (eye(nF) - ones(nF)*b/(a+b*nF))/a;
                            
                            [R(1,2),R(2,1)] = deal(hha(1)+hha(2));
                            det_Va          = det_Va*det(R)^n_nFMZ;
                            invMZR{nF}      = sparse(pinv(R));
                            
                            %%% CE model
                            R0              = eye(nF)*hh0(3) + ones(nF)*hh0(2);
                            det_V0          = det_V0*det(R0)^n_nF;
                            invR0{nF}       = (eye(nF) - ones(nF)*hh0(2)/(hh0(3)+hh0(2)*nF))/hh0(3);
                        end
                        
                        % ACE model
                        Ua = blkdiag(invOTR{nFam});
                        for nF = Perm_label(1:nMZF)
                            IndF = nFam(nF);
                            i0   = sum(nFam(1:nF-1));
                            Ua(i0+[1:IndF], i0+[1:IndF]) = invMZR{IndF};
                        end
                        Za  = pinv(X'*Ua*X);
                        Pa  = Ua-Ua*X*Za*X'*Ua;
                        rla = -(iY'*Pa*iY/sig2+log(det_Va)+log(det(X'*Ua*X)))/2;
                        
                        % CE model
                        U0  = blkdiag(invR0{nFam});
                        Z0  = pinv(X'*U0*X);
                        P0  = U0-U0*X*Z0*X'*U0;
                        rl0 = -(iY'*P0*iY/sig2+log(det_V0)+log(det(X'*U0*X)))/2;
                        
                        % Compute LRT statistic
                        Stats(i) = max(0,-2*(rl0-rla));
                    end
                    
                end
                
            end
            
            % Compute maximum LRT statistic
            mT = max(Stats);
            
            if length(Vs.Dim)==1
                Dim = [Vs.Dim 1];
            else
                Dim = Vs.Dim;
            end
            
            Stats  = reshape(Stats,Dim);
            Tstats = reshape(ACEfit_Par.Stats,Dim);
            
            uCount = zeros(Dim);
            uCount(Stats>=Tstats) = 1;
            
            %
            % Cluster inference
            %
            if ACEfit_Par.Vs.ClustInf
                
                CFT = spm_invXcdf(1-2*ACEfit_Par.alpha_CFT,1);
                
                % Apply cluster-forming threshold (CFT)
                Stats(Stats<CFT) = 0;
                
                % Cluster-forming procedure
                [L,NUM] = spm_bwlabel(Stats,18);
                
                if NUM>0
                    cluster_size = zeros(NUM,1);
                    cluster_mass = zeros(NUM,1);
                    for i=1:NUM
                        cluster_size(i) = length(find(L(:)==i));
                        cluster_mass(i) = sum(Stats(L(:)==i));
                    end
                    % Compute maximum cluster statistics
                    mK = max(cluster_size);
                    mM = max(cluster_mass);
                end
                
            end
            
        end
        
    case 'AE'
        
        hh           = zeros(nElm,3);
        hh(I_data,:) = hhAE0./repmat(sum(hhAE0,2),1,3);
        
        %
        % Image-wise inference
        %
        if ~ACEfit_Par.NoImg
            
            %
            % Likelihood ratio test
            %
            
            % Initialization
            Stats = zeros(nElm,1);
            
            P0    = eye(n)-X*pinv(X);
            
            for i = I_data
                
                iY   = Y(i,:)';
                hha  = hh(i,:)';
                sig2 = Sigma2(i);
                
                if (hha(3)<1e-5)
                    % Zero environmental variance... give up
                    BadElm   = BadElm + 1;
                    Stats(i) = 0;
                else
                    
                    % Parameter A is non-zero... compute log-likelihood...
                    index = ACEcode(hha);
                    if ( index==ACEcode([1 0 1]) )
                        
                        det_Va                = 1;
                        [invMZR,invOTR]       = deal(cell(max(nFam),1));
                        [invMZR{1},invOTR{1}] = deal(1);
                        
                        for nF = 2:max(nFam)
                            n_nF   = sum(nFam==nF);
                            n_nFMZ = sum(nFam(Perm_label(1:nMZF))==nF);
                            
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
                        Ua = blkdiag(invOTR{nFam});
                        for nF = Perm_label(1:nMZF)
                            IndF = nFam(nF);
                            i0   = sum(nFam(1:nF-1));
                            Ua(i0+[1:IndF], i0+[1:IndF]) = invMZR{IndF};
                        end
                        Za  = pinv(X'*Ua*X);
                        Pa  = Ua-Ua*X*Za*X'*Ua;
                        rla = -(iY'*Pa*iY/sig2+log(det_Va)+log(det(X'*Ua*X)))/2;
                        
                        % E model
                        rl0 = -(iY'*P0*iY/sig2+log(det(X'*X)))/2;
                        
                        % Compute LRT statistic
                        Stats(i) = max(0,-2*(rl0-rla));
                    end
                    
                end
                
            end
            
            % Compute maximum LRT statistic
            mT = max(Stats);
            
            if length(Vs.Dim)==1
                Dim = [Vs.Dim 1];
            else
                Dim = Vs.Dim;
            end
            
            Stats  = reshape(Stats,Dim);
            Tstats = reshape(ACEfit_Par.Stats,Dim);
            
            uCount = zeros(Dim);
            uCount(Stats>=Tstats) = 1;
            
            %
            % Cluster inference
            %
            if ACEfit_Par.Vs.ClustInf
                
                CFT = spm_invXcdf(1-2*ACEfit_Par.alpha_CFT,1);
                
                % Apply cluster-forming threshold (CFT)
                Stats(Stats<CFT) = 0;
                
                % Cluster-forming procedure
                [L,NUM] = spm_bwlabel(Stats,18);
                
                if NUM>0
                    cluster_size = zeros(NUM,1);
                    cluster_mass = zeros(NUM,1);
                    for i=1:NUM
                        cluster_size(i) = length(find(L(:)==i));
                        cluster_mass(i) = sum(Stats(L(:)==i));
                    end
                    % Compute maximum cluster statistics
                    mK = max(cluster_size);
                    mM = max(cluster_mass);
                end
                
            end
            
        end
        
end

% Summary statistics
h2       = hh(I_data,1);
Sigma2   = Sigma2(I_data);

Aperm    = Sigma2.*h2;
wh2      = mean(Aperm)/mean(Sigma2);

quartile = quantile(h2,[0.50 0.75]);
medh2    = quartile(1);
q3h2     = quartile(2);
meanh2   = mean(h2);
mGmedh2  = mean(h2(h2>medh2));
mGq3h2   = mean(h2(h2>q3h2));

% Summary statistics: meanh2, wh2, median, q3, mean(h2>median), mean(h2>q3)
SummaryA = [meanh2; wh2; medh2; q3h2; mGmedh2; mGq3h2];

return



