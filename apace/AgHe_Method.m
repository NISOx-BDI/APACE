function AgHe_Method(ACEfit_Par,varargin)
%
% Permutation & bootstrap inferences for AgHe (aka Steve's) method
%
% AgHe_Method(ACEfit_Par)                        - Save results
% AgHe_Method(ACEfit_Par,'_Norm')                - Save results with suffix
%                                                  for permutation & bootstrap
% AgHe_Method(ACEfit_Par,'_Norm',Palpha,Balpha)  - Save results with the suffix,
%                                                  also set Perm & Bootstrap alpha
%_______________________________________________________________________
% Version: http://github.com/nicholst/APACE/tree/$Format:%h$
%          $Format:%ci$

if size(ACEfit_Par.Y,1)==1
    error('Aggregate heritability is designed for multiple phenotypes!')
end

if nargin<2
    ResSuf = '';
else
    ResSuf = varargin{1};
end
if nargin<3 || isempty(varargin{2})
    Palpha = 0.05;
else
    Palpha = varargin{2};
end
if nargin<4 || isempty(varargin{3})
    Balpha = 0.05;
else
    Balpha = varargin{3};
end

nMZF  = ACEfit_Par.nMZF;
nDZF  = ACEfit_Par.nDZF;
nSibF = ACEfit_Par.nSibF;

nFam  = ACEfit_Par.nFam;
cFam  = [1; cumsum(nFam)+1];
cFam  = cFam(1:end-1);

%
% Compute the OLS residuals, with appropriate normalisation
%
Res = ACEfit_Par.Y' - ACEfit_Par.X*pinv(ACEfit_Par.X)*ACEfit_Par.Y';

if ( ~isfield(ACEfit_Par,'AggNlz') || isempty(ACEfit_Par.AggNlz) )
    ACEfit_Par.AggNlz = 0;
end

switch(ACEfit_Par.AggNlz)
    case 1
        for j=1:size(Res,2)
            Res(:,j) = Res(:,j)/std(Res(:,j));
        end
    case 2
        for j=1:size(Res,2)
            Res(:,j) = Res(:,j) + ACEfit_Par.Ymean(j);
        end
    case 3
        for j=1:size(Res,2)
            Res(:,j) = (Res(:,j)+ACEfit_Par.Ymean(j))/std(Res(:,j));
        end
end

% Apply mask
Res = Res(:,ACEfit_Par.I_data);


%
% Select sibling pairs from non-DZ families
%
rng('shuffle');

% Select sibling pairs from MZ families
MZFindex = zeros(nMZF,2);

iMZF2 = nFam(1:nMZF)>2;
nMZF2 = sum(iMZF2);

label_MZ = reshape(repmat(1:2,nMZF2,1)',2*nMZF2,1);
MZindex  = reshape(label_MZ(randperm(2*nMZF2)),2,nMZF2)';
MZFindex(iMZF2,1) = MZindex(:,1);

MZFindex(nFam(1:nMZF)==3,2) = 3;

iMZF3 = nFam(1:nMZF)>3;
nMZF3 = sum(iMZF3);

Isibs  = nFam(iMZF3)-2;
MZFsib = zeros(nMZF3,1);
for k=1:nMZF3
    MZFsib(k) = randperm(Isibs(k),1)+2;
end
MZFindex(iMZF3,2) = MZFsib;

MZFindex(iMZF2,:) = MZFindex(iMZF2,:)+repmat(cFam(iMZF2),1,2)-1;
MZFindex          = MZFindex(iMZF2,:);

% Select sibling pairs from non-twin families
SibFindex = zeros(nSibF,2);
nFamSibF  = nFam(nMZF+nDZF+[1:nSibF]);

SibFindex(nFamSibF==2,1) = 1;
SibFindex(nFamSibF==2,2) = 2;

iSibF2 = nFamSibF>2;
nSibF2 = sum(iSibF2);

Isibs   = nFamSibF(iSibF2);
SibFsib = zeros(nSibF2,2);
for k=1:nSibF2
    SibFsib(k,:) = randperm(Isibs(k),2);
end
SibFindex(iSibF2,:) = SibFsib;

SibFindex = SibFindex+repmat(cFam(nMZF+nDZF+[1:nSibF]),1,2)-1;

% Index of sibling pairs from non-DZ families
SibPindex = [MZFindex; SibFindex];
% Number of sibling pairs
nSibP = size(SibPindex,1);

% Generate permutation and bootstrapping labels for rDZ & rSib
if nSibP>0
    if ACEfit_Par.nPerm>0
        [Perm_DZSib] = CreatePerm(nDZF,nSibP,ACEfit_Par.nPerm);
    end
    if ACEfit_Par.nBoot>0
        [Boot_DZSib] = CreateBoot(nDZF,nSibP,0,ACEfit_Par.nBoot);
    end
end


%
% Compute rMZo, rDZo, rSibo for the original data
%
ResY1 = [Res(cFam(1:nMZF+nDZF),:);   Res(SibPindex(:,1),:)];
ResY2 = [Res(cFam(1:nMZF+nDZF)+1,:); Res(SibPindex(:,2),:)];

rMZo = zeros(nMZF,1);
for j=1:nMZF
    MZ1     = ResY1(j,:);
    MZ2     = ResY2(j,:);
    rMZ     = corrcoef(MZ1,MZ2);
    rMZo(j) = rMZ(1,2);
end
rDZo = zeros(nDZF,1);
for j=1:nDZF
    DZ1     = ResY1(nMZF+j,:);
    DZ2     = ResY2(nMZF+j,:);
    rDZ     = corrcoef(DZ1,DZ2);
    rDZo(j) = rDZ(1,2);
end
% Point estimates of rMZ-rDZ
Est_MZDZ = mean(rMZo) - mean(rDZo);

if nSibP>0
    rSibo = zeros(nSibP,1);
    for j=1:nSibP
        Sib1     = ResY1(nMZF+nDZF+j,:);
        Sib2     = ResY2(nMZF+nDZF+j,:);
        rSib     = corrcoef(Sib1,Sib2);
        rSibo(j) = rSib(1,2);
    end
    % Point estimates of rDZ-rSib
    Est_DZSib = mean(rDZo) - mean(rSibo);
end


f = get(0,'Children');
if isempty(f)
    f = 1;
else
    f = length(f)+1;
    % f = max(f)+1;
end

SetFig(f);
% AgHe box plots
if nSibP>0
    Z = [rMZo;                  rDZo;                 rSibo                  ];
    g = [zeros(length(rMZo),1); ones(length(rDZo),1); 2*ones(length(rSibo),1)];
    boxplot(Z,g);
    title(sprintf('Boxplots of Correlations of Aggregate Heritabilities'));
    set(gca, 'XTick', [1; 2; 3]);
    set(gca, 'XTickLabel', {'MZ'; 'DZ'; 'SIB'});
    print('-dpdf',fullfile(ACEfit_Par.ResDir,'AgHe_CorrPlots.pdf'));
    
    save(fullfile(ACEfit_Par.ResDir,['Ests_AgHe' ResSuf]),...
                                     'Est_MZDZ','Est_DZSib',...
                                     'rMZo', 'rDZo', 'rSibo');
else
    Z = [rMZo;                  rDZo;               ];
    g = [zeros(length(rMZo),1); ones(length(rDZo),1)];
    boxplot(Z,g);
    title(sprintf('Boxplots of Correlations of Aggregate Heritabilities'));
    set(gca, 'XTick', [1; 2]);
    set(gca, 'XTickLabel', {'MZ'; 'DZ'});
    print('-dpdf',fullfile(ACEfit_Par.ResDir,'AgHe_CorrPlots.pdf'));
    
    save(fullfile(ACEfit_Par.ResDir,['Ests_AgHe' ResSuf]),...
                                     'Est_MZDZ',...
                                     'rMZo', 'rDZo');
end


%
% Permutations
%

if ACEfit_Par.nPerm>0
    
    N = ACEfit_Par.nPerm+1;
    
    corrMZ1 = zeros(N,nMZF);
    corrDZ1 = zeros(N,nDZF);
    
    corrMZ1(end,:) = rMZo;
    corrDZ1(end,:) = rDZo;
    
    % Compute rMZ & rDZ for all permutations
    for i=1:ACEfit_Par.nPerm
        Perm_label    = ACEfit_Par.Perm_index(i,:);
        corrMZDZ      = [corrMZ1(end,:) corrDZ1(end,:)];
        corrMZDZ      = corrMZDZ(Perm_label);
        corrMZ1(i,:)  = corrMZDZ(1:nMZF);
        corrDZ1(i,:)  = corrMZDZ(nMZF+[1:nDZF]);
    end
    
    % 1st MZ vs DZ statistic: 2-sample t-statistic
    t_stats_MZ_DZ  = (mean(corrMZ1,2)-mean(corrDZ1,2))./sqrt(var(corrMZ1,0,2)/nMZF+var(corrDZ1,0,2)/nDZF);
    
    % 2nd MZ vs DZ statistic: mean difference
    CorrDiff_MZ_DZ = mean(corrMZ1,2)-mean(corrDZ1,2);
    

    if nSibP>0
        
        corrDZ2  = zeros(N,nDZF);
        corrSib2 = zeros(N,nSibP);
        
        corrDZ2(end,:)  = rDZo;
        corrSib2(end,:) = rSibo;
        
        % Compute rDZ & rSib for all permutations
        for i=1:ACEfit_Par.nPerm
            Perm_label    = Perm_DZSib(i,:);
            corrDZSib     = [corrDZ2(end,:) corrSib2(end,:)];
            corrDZSib     = corrDZSib(Perm_label);
            corrDZ2(i,:)  = corrDZSib(1:nDZF);
            corrSib2(i,:) = corrDZSib(nDZF+[1:nSibP]);
        end
        
        % 1st DZ vs Sib statistic: 2-sample t-statistic
        t_stats_DZ_Sib  = (mean(corrDZ2,2)-mean(corrSib2,2))./sqrt(var(corrDZ2,0,2)/nDZF+var(corrSib2,0,2)/nSibP);
        
        % 2nd DZ vs Sib statistic: mean difference
        CorrDiff_DZ_Sib = mean(corrDZ2,2)-mean(corrSib2,2);
        
    end
    
    
    %
    % Plot the empirical distributions of the statistics
    %
    ctl_val_index = N-floor(N*Palpha);
    
    %
    % 1st MZ vs DZ statistic: 2-sample t-statistic
    %
    f = f+1;
    SetFig(f);
    [st_stats_MZ_DZ,~]    = sort(t_stats_MZ_DZ);
    ctl_val_t_stats_MZ_DZ = st_stats_MZ_DZ(ctl_val_index);
    [f1,x1]               = hist(st_stats_MZ_DZ,100);
    n_t_stats_MZ_DZ       = min(find(st_stats_MZ_DZ(:)==t_stats_MZ_DZ(N)));
    bar(x1,f1/N);
    
    % Calculate the permutation-based p-value
    p_t_stats_MZ_DZ = (N-n_t_stats_MZ_DZ+1)/N;
    
    title(sprintf('H0 dist of two-sample t-statistic for rMZ & rDZ, P-value=%.3f',p_t_stats_MZ_DZ));
    xlabel('t-statistic for rMZ vs. rDZ');
    ylabel('frequency');
    yLimits = get(gca,'YLim');
    line([t_stats_MZ_DZ(N) t_stats_MZ_DZ(N)],[0 yLimits(2)],'Marker','.','Color','green');
    line([ctl_val_t_stats_MZ_DZ ctl_val_t_stats_MZ_DZ],[0 yLimits(2)],'Marker','.','LineStyle','-.','Color','red');
    text(t_stats_MZ_DZ(N),0.8*yLimits(2),num2str(t_stats_MZ_DZ(N)))
    text(ctl_val_t_stats_MZ_DZ,0.9*yLimits(2),num2str(ctl_val_t_stats_MZ_DZ))
    print('-dpdf',fullfile(ACEfit_Par.ResDir,['H0dist_Tstat_rMZ_rDZ' ResSuf '.pdf']));
    
    %
    % 2nd MZ vs DZ statistic: mean difference
    %
    f = f+1;
    SetFig(f);
    [sCorrDiff_MZ_DZ,~]    = sort(CorrDiff_MZ_DZ);
    ctl_val_CorrDiff_MZ_DZ = sCorrDiff_MZ_DZ(ctl_val_index);
    [f2,x2]                = hist(sCorrDiff_MZ_DZ,100);
    n_CorrDiff_MZ_DZ       = min(find(sCorrDiff_MZ_DZ(:)==CorrDiff_MZ_DZ(N)));
    bar(x2,f2/N);
    
    % Calculate the permutation-based p-value
    p_CorrDiff_MZ_DZ = (N-n_CorrDiff_MZ_DZ+1)/N;
    
    title(sprintf('H0 dist of mean difference between rMZ & rDZ, P-value=%.3f',p_CorrDiff_MZ_DZ));
    xlabel('E(rMZ) - E(rDZ)');
    ylabel('frequency');
    yLimits = get(gca,'YLim');
    line([CorrDiff_MZ_DZ(N) CorrDiff_MZ_DZ(N)],[0 yLimits(2)],'Marker','.','Color','green');
    line([ctl_val_CorrDiff_MZ_DZ ctl_val_CorrDiff_MZ_DZ],[0 yLimits(2)],'Marker','.','LineStyle','-.','Color','red');
    text(CorrDiff_MZ_DZ(N),0.8*yLimits(2),num2str(CorrDiff_MZ_DZ(N)))
    text(ctl_val_CorrDiff_MZ_DZ,0.9*yLimits(2),num2str(ctl_val_CorrDiff_MZ_DZ))
    print('-dpdf',fullfile(ACEfit_Par.ResDir,['H0dist_Diff_rMZ_rDZ' ResSuf '.pdf']));
    
    Pvals_MZDZ  = [p_t_stats_MZ_DZ  p_CorrDiff_MZ_DZ]';
    
    
    if nSibP>0
        
        %
        % 1st DZ vs Sib statistic: 2-sample t-statistic
        %
        f = f+1;
        SetFig(f);
        [st_stats_DZ_Sib,~]    = sort(t_stats_DZ_Sib);
        ctl_val_t_stats_DZ_Sib = st_stats_DZ_Sib(ctl_val_index);
        [f1,x1]                = hist(st_stats_DZ_Sib,100);
        n_t_stats_DZ_Sib       = min(find(st_stats_DZ_Sib(:)==t_stats_DZ_Sib(N)));
        bar(x1,f1/N);
        
        % Calculate the permutation-based p-value
        p_t_stats_DZ_Sib = (N-n_t_stats_DZ_Sib+1)/N;
        
        title(sprintf('H0 dist of two-sample t-statistic for rDZ & rSib, P-value=%.3f',p_t_stats_DZ_Sib));
        xlabel('t-statistic for rDZ vs. rSib');
        ylabel('frequency');
        yLimits = get(gca,'YLim');
        line([t_stats_DZ_Sib(N) t_stats_DZ_Sib(N)],[0 yLimits(2)],'Marker','.','Color','green');
        line([ctl_val_t_stats_DZ_Sib ctl_val_t_stats_DZ_Sib],[0 yLimits(2)],'Marker','.','LineStyle','-.','Color','red');
        text(t_stats_DZ_Sib(N),0.8*yLimits(2),num2str(t_stats_DZ_Sib(N)))
        text(ctl_val_t_stats_DZ_Sib,0.9*yLimits(2),num2str(ctl_val_t_stats_DZ_Sib))
        print('-dpdf',fullfile(ACEfit_Par.ResDir,['H0dist_Tstat_rDZ_rSib' ResSuf '.pdf']));
        
        %
        % 2nd DZ vs Sib statistic: mean difference
        %
        f = f+1;
        SetFig(f);
        [sCorrDiff_DZ_Sib,~]    = sort(CorrDiff_DZ_Sib);
        ctl_val_CorrDiff_DZ_Sib = sCorrDiff_DZ_Sib(ctl_val_index);
        [f2,x2]                 = hist(sCorrDiff_DZ_Sib,100);
        n_CorrDiff_DZ_Sib       = min(find(sCorrDiff_DZ_Sib(:)==CorrDiff_DZ_Sib(N)));
        bar(x2,f2/N);
        
        % Calculate the permutation-based p-value
        p_CorrDiff_DZ_Sib = (N-n_CorrDiff_DZ_Sib+1)/N;
        
        title(sprintf('H0 dist of mean difference between rDZ & rSib, P-value=%.3f',p_CorrDiff_DZ_Sib));
        xlabel('E(rDZ) - E(rSib)');
        ylabel('frequency');
        yLimits = get(gca,'YLim');
        line([CorrDiff_DZ_Sib(N) CorrDiff_DZ_Sib(N)],[0 yLimits(2)],'Marker','.','Color','green');
        line([ctl_val_CorrDiff_DZ_Sib ctl_val_CorrDiff_DZ_Sib],[0 yLimits(2)],'Marker','.','LineStyle','-.','Color','red');
        text(CorrDiff_DZ_Sib(N),0.8*yLimits(2),num2str(CorrDiff_DZ_Sib(N)))
        text(ctl_val_CorrDiff_DZ_Sib,0.9*yLimits(2),num2str(ctl_val_CorrDiff_DZ_Sib))
        print('-dpdf',fullfile(ACEfit_Par.ResDir,['H0dist_Diff_rDZ_rSib' ResSuf '.pdf']));
        
        Pvals_DZSib = [p_t_stats_DZ_Sib p_CorrDiff_DZ_Sib]';
        
        % Save p-values
        save(fullfile(ACEfit_Par.ResDir,['Pvals_AgHe' ResSuf]),'Pvals_MZDZ','Pvals_DZSib');
    else
        save(fullfile(ACEfit_Par.ResDir,['Pvals_AgHe' ResSuf]),'Pvals_MZDZ');
    end
    
end


%
% Bootstrapping
%

if ACEfit_Par.nBoot>0
    
    nB = ACEfit_Par.nBoot+1;
    
    Boot_rMZ1  = zeros(nB,nMZF);
    Boot_rDZ1  = zeros(nB,nDZF);
    
    Boot_rMZ1(end,:)  = rMZo;
    Boot_rDZ1(end,:)  = rDZo;
    
    % Compute rMZ & rDZ for all bootstrap replicates
    Boot_MZDZ = ACEfit_Par.Boot_index(:,1:nMZF+nDZF);
    for i=1:ACEfit_Par.nBoot
        Boot_label      = Boot_MZDZ(i,:);
        Boot_rMZ1(i,:)  = Boot_rMZ1(end,Boot_label(1:nMZF));
        Boot_rDZ1(i,:)  = Boot_rDZ1(end,Boot_label(nMZF+[1:nDZF])-nMZF);
    end
    
    % Bootstrapping MZ vs DZ mean differences
    Boot_diff_MZ_DZ  = mean(Boot_rMZ1,2)-mean(Boot_rDZ1,2);
    
    % Confidence level 1-Balpha
    CI_MZDZ  = zeros(1,2);
    [CI_MZDZ(1), CI_MZDZ(2)] = Bootstrap_CI(Boot_diff_MZ_DZ, Balpha);
    
    if nSibP>0
        Boot_rDZ2  = zeros(nB,nDZF);
        Boot_rSib2 = zeros(nB,nSibP);
        
        Boot_rDZ2(end,:)  = rDZo;
        Boot_rSib2(end,:) = rSibo;
        
        % Compute rDZ & rSib for all bootstrap replicates
        for i=1:ACEfit_Par.nBoot
            Boot_label      = Boot_DZSib(i,:);
            Boot_rDZ2(i,:)  = Boot_rDZ2(end,Boot_label(1:nDZF));
            Boot_rSib2(i,:) = Boot_rSib2(end,Boot_label(nDZF+[1:nSibP])-nDZF);
        end
        
        % Bootstrapping DZ vs Sib mean differences
        Boot_diff_DZ_Sib = mean(Boot_rDZ2,2)-mean(Boot_rSib2,2);
        
        % Confidence level 1-Balpha
        CI_DZSib = zeros(1,2);
        [CI_DZSib(1), CI_DZSib(2)] = Bootstrap_CI(Boot_diff_DZ_Sib, Balpha);
        
        % Save CI's
        save(fullfile(ACEfit_Par.ResDir,['CIs_AgHe' ResSuf]),'CI_MZDZ','CI_DZSib');
        
    else
        save(fullfile(ACEfit_Par.ResDir,['CIs_AgHe' ResSuf]),'CI_MZDZ');
    end

end


return

function SetFig(f)

figure(f)
set(gcf,'PaperPosition',[0 0 10 6])
set(gcf,'PaperSize',[10 6])

return
