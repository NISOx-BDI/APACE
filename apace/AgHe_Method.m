function AgHe_Method(ACEfit_Par,Palpha,Balpha,varargin)
%
% Permutation & bootstrap inferences for AgHe (aka Steve's) method
%
% AgHe_Method(ACEfit_Par,Palpha,Balpha)          - Saves results
% AgHe_Method(ACEfit_Par,Palpha,Balpha,'_Norm')  - Saves results with the suffix
%

if nargin<=3
    ResSuf = '';
else
    % Results suffix
    ResSuf = varargin{1};
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
% if ACEfit_Par.AggNlz
%     Res = ACEfit_Par.Y' - ACEfit_Par.X*pinv(ACEfit_Par.X)*ACEfit_Par.Y';
% else
%     Res = ACEfit_Par.Y';
% end

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
if ( nDZF==0 || nSibP==0 )
    error('Cannot make permutation inference as there are no DZ twins or sibling pairs!!')
end
if ACEfit_Par.nPerm>0
    [Perm_DZSib] = CreatePerm(nDZF,nSibP,ACEfit_Par.nPerm);
end
if ACEfit_Par.nBoot>0
    [Boot_DZSib] = CreateBoot(nDZF,nSibP,0,ACEfit_Par.nBoot);
end


%
% Compute rMZo, rDZo, rSibo for the original data
%
ResY1 = [Res(cFam(1:nMZF+nDZF),:);   Res(SibPindex(:,1),:)];
ResY2 = [Res(cFam(1:nMZF+nDZF)+1,:); Res(SibPindex(:,2),:)];

rMZo  = zeros(nMZF,1);
rDZo  = zeros(nDZF,1);
rSibo = zeros(nSibP,1);

for j=1:nMZF
    MZ1      = ResY1(j,:);
    MZ2      = ResY2(j,:);
    rMZ      = corrcoef(MZ1,MZ2);
    rMZo(j)  = rMZ(1,2);
end
for j=1:nDZF
    DZ1      = ResY1(nMZF+j,:);
    DZ2      = ResY2(nMZF+j,:);
    rDZ      = corrcoef(DZ1,DZ2);
    rDZo(j)  = rDZ(1,2);
end
for j=1:nSibP
    Sib1     = ResY1(nMZF+nDZF+j,:);
    Sib2     = ResY2(nMZF+nDZF+j,:);
    rSib     = corrcoef(Sib1,Sib2);
    rSibo(j) = rSib(1,2);
end

% Point estimates of rMZ-rDZ & rDZ-rSib
Est_MZDZ  = mean(rMZo) - mean(rDZo);
Est_DZSib = mean(rDZo) - mean(rSibo);


f=max(get(0,'Children'))+1;
if isempty(f); f=1; end
figure(f); 
% AgHe box plots
Z = [rMZo;                  rDZo;                 rSibo                  ];
g = [zeros(length(rMZo),1); ones(length(rDZo),1); 2*ones(length(rSibo),1)];
boxplot(Z,g);
title(sprintf('Boxplots of Correlations of Aggregate Heritabilities'));
set(gca, 'XTick', [1; 2; 3]);
set(gca, 'XTickLabel', {'MZ'; 'DZ'; 'SIB'});
set(gcf,'PaperPosition',[0 0 10 6])
set(gcf,'PaperSize',[10 6])
print('-dpdf',fullfile(ACEfit_Par.ResDir, 'AgHe_CorrPlots.pdf'));




save(fullfile(ACEfit_Par.ResDir,['Ests_AgHe' ResSuf]),...
     'Est_MZDZ','Est_DZSib',...
     'rMZo', 'rDZo', 'rSibo');


%
% Permutations
%

if ACEfit_Par.nPerm>0
    
    N = ACEfit_Par.nPerm+1;
    
    corrMZ1  = zeros(N,nMZF);
    corrDZ1  = zeros(N,nDZF);
    
    corrDZ2  = zeros(N,nDZF);
    corrSib2 = zeros(N,nSibP);
    
    corrMZ1(end,:)  = rMZo;
    corrDZ1(end,:)  = rDZo;
    
    corrDZ2(end,:)  = rDZo;
    corrSib2(end,:) = rSibo;
    
    % Compute rMZ & rDZ for all permutations
    for i=1:ACEfit_Par.nPerm
        Perm_label    = ACEfit_Par.Perm_index(i,:);
        corrMZDZ      = [corrMZ1(end,:) corrDZ1(end,:)];
        corrMZDZ      = corrMZDZ(Perm_label);
        corrMZ1(i,:)  = corrMZDZ(1:nMZF);
        corrDZ1(i,:)  = corrMZDZ(nMZF+[1:nDZF]);
    end
    
    % Compute rDZ & rSib for all permutations
    for i=1:ACEfit_Par.nPerm
        Perm_label    = Perm_DZSib(i,:);
        corrDZSib     = [corrDZ2(end,:) corrSib2(end,:)];
        corrDZSib     = corrDZSib(Perm_label);
        corrDZ2(i,:)  = corrDZSib(1:nDZF);
        corrSib2(i,:) = corrDZSib(nDZF+[1:nSibP]);
    end
    
    % 1st statistic: 2-sample t-statistic
    t_stats_MZ_DZ   = (mean(corrMZ1,2)-mean(corrDZ1,2))./sqrt(var(corrMZ1,0,2)/nMZF+var(corrDZ1,0,2)/nDZF);
    t_stats_DZ_Sib  = (mean(corrDZ2,2)-mean(corrSib2,2))./sqrt(var(corrDZ2,0,2)/nDZF+var(corrSib2,0,2)/nSibP);
    
    % 2nd statistic: mean difference
    CorrDiff_MZ_DZ  = mean(corrMZ1,2)-mean(corrDZ1,2);
    CorrDiff_DZ_Sib = mean(corrDZ2,2)-mean(corrSib2,2);
    
    
    %
    % Plot the empirical distributions of the statistics
    %
    ctl_val_index = N-floor(N*Palpha);
    
    f = 0;
    
    %
    % 1st statistic: 2-sample t-statistic
    %
    f = f+1;
    figure(f);
    [st_stats_MZ_DZ,~]    = sort(t_stats_MZ_DZ);
    ctl_val_t_stats_MZ_DZ = st_stats_MZ_DZ(ctl_val_index);
    [f1,x1]               = hist(st_stats_MZ_DZ,100);
    n_t_stats_MZ_DZ       = min(find(st_stats_MZ_DZ(:)==t_stats_MZ_DZ(N)));
    bar(x1,f1/N);
    
    % Calculate the permutation-based p-value
    p_t_stats_MZ_DZ = (N-n_t_stats_MZ_DZ+1)/N;
    
    set(gca, 'xtick', 0)
    if ctl_val_t_stats_MZ_DZ==t_stats_MZ_DZ(N)
        set(gca, 'xtick', ctl_val_t_stats_MZ_DZ)
    else
        set(gca, 'xtick', sort([ctl_val_t_stats_MZ_DZ t_stats_MZ_DZ(N)]))
    end
    
    title(sprintf('Empirical distribution of two-sample t-statistic for rMZ & rDZ, P-value=%.3f',p_t_stats_MZ_DZ));
    xlabel('t-statistic for rMZ vs. rDZ');
    ylabel('frequency');
    hold on;
    yLimits = get(gca,'YLim');
    line([t_stats_MZ_DZ(N) t_stats_MZ_DZ(N)],[0 yLimits(2)],'Marker','.','Color','green');
    line([ctl_val_t_stats_MZ_DZ ctl_val_t_stats_MZ_DZ],[0 yLimits(2)],'Marker','.','LineStyle','-.','Color','red');
    hold off;
    set(gcf,'PaperPosition',[0 0 10 6])
    set(gcf,'PaperSize',[10 6])
    print('-dpdf',fullfile(ACEfit_Par.ResDir,['Tstat_rMZ_rDZ' ResSuf '.pdf']));
    
    
    f = f+1;
    figure(f);
    [st_stats_DZ_Sib,~]    = sort(t_stats_DZ_Sib);
    ctl_val_t_stats_DZ_Sib = st_stats_DZ_Sib(ctl_val_index);
    [f1,x1]                = hist(st_stats_DZ_Sib,100);
    n_t_stats_DZ_Sib       = min(find(st_stats_DZ_Sib(:)==t_stats_DZ_Sib(N)));
    bar(x1,f1/N);
    
    % Calculate the permutation-based p-value
    p_t_stats_DZ_Sib = (N-n_t_stats_DZ_Sib+1)/N;
    
    set(gca, 'xtick', 0)
    if ctl_val_t_stats_DZ_Sib==t_stats_DZ_Sib(N)
        set(gca, 'xtick', ctl_val_t_stats_DZ_Sib)
    else
        set(gca, 'xtick', sort([ctl_val_t_stats_DZ_Sib t_stats_DZ_Sib(N)]))
    end
    
    title(sprintf('Empirical distribution of two-sample t-statistic for rDZ & rSib, P-value=%.3f',p_t_stats_DZ_Sib));
    xlabel('t-statistic for rDZ vs. rSib');
    ylabel('frequency');
    hold on;
    yLimits = get(gca,'YLim');
    line([t_stats_DZ_Sib(N) t_stats_DZ_Sib(N)],[0 yLimits(2)],'Marker','.','Color','green');
    line([ctl_val_t_stats_DZ_Sib ctl_val_t_stats_DZ_Sib],[0 yLimits(2)],'Marker','.','LineStyle','-.','Color','red');
    hold off;
    set(gcf,'PaperPosition',[0 0 10 6])
    set(gcf,'PaperSize',[10 6])
    print('-dpdf',fullfile(ACEfit_Par.ResDir,['Tstat_rDZ_rSib' ResSuf '.pdf']));
    
    
    %
    % 2nd statistic: mean difference
    %
    f = f+1;
    figure(f);
    [sCorrDiff_MZ_DZ,~]    = sort(CorrDiff_MZ_DZ);
    ctl_val_CorrDiff_MZ_DZ = sCorrDiff_MZ_DZ(ctl_val_index);
    [f2,x2]                = hist(sCorrDiff_MZ_DZ,100);
    n_CorrDiff_MZ_DZ       = min(find(sCorrDiff_MZ_DZ(:)==CorrDiff_MZ_DZ(N)));
    bar(x2,f2/N);
    
    % Calculate the permutation-based p-value
    p_CorrDiff_MZ_DZ = (N-n_CorrDiff_MZ_DZ+1)/N;
    
    set(gca, 'xtick', 0)
    if ctl_val_CorrDiff_MZ_DZ==CorrDiff_MZ_DZ(N)
        set(gca, 'xtick', ctl_val_CorrDiff_MZ_DZ)
    else
        set(gca, 'xtick', sort([ctl_val_CorrDiff_MZ_DZ CorrDiff_MZ_DZ(N)]))
    end
    
    title(sprintf('Empirical distribution of mean difference between rMZ & rDZ, P-value=%.3f',p_CorrDiff_MZ_DZ));
    xlabel('E(rMZ) - E(rDZ)');
    ylabel('frequency');
    hold on;
    yLimits = get(gca,'YLim');
    line([CorrDiff_MZ_DZ(N) CorrDiff_MZ_DZ(N)],[0 yLimits(2)],'Marker','.','Color','green');
    line([ctl_val_CorrDiff_MZ_DZ ctl_val_CorrDiff_MZ_DZ],[0 yLimits(2)],'Marker','.','LineStyle','-.','Color','red');
    hold off;
    set(gcf,'PaperPosition',[0 0 10 6])
    set(gcf,'PaperSize',[10 6])
    print('-dpdf',fullfile(ACEfit_Par.ResDir,['Diff_rMZ_rDZ' ResSuf '.pdf']));
    
    
    f = f+1;
    figure(f);
    [sCorrDiff_DZ_Sib,~]    = sort(CorrDiff_DZ_Sib);
    ctl_val_CorrDiff_DZ_Sib = sCorrDiff_DZ_Sib(ctl_val_index);
    [f2,x2]                 = hist(sCorrDiff_DZ_Sib,100);
    n_CorrDiff_DZ_Sib       = min(find(sCorrDiff_DZ_Sib(:)==CorrDiff_DZ_Sib(N)));
    bar(x2,f2/N);
    
    % Calculate the permutation-based p-value
    p_CorrDiff_DZ_Sib = (N-n_CorrDiff_DZ_Sib+1)/N;
    
    set(gca, 'xtick', 0)
    if ctl_val_CorrDiff_DZ_Sib==CorrDiff_DZ_Sib(N)
        set(gca, 'xtick', ctl_val_CorrDiff_DZ_Sib)
    else
        set(gca, 'xtick', sort([ctl_val_CorrDiff_DZ_Sib CorrDiff_DZ_Sib(N)]))
    end
    
    title(sprintf('Empirical distribution of mean difference between rDZ & rSib, P-value=%.3f',p_CorrDiff_DZ_Sib));
    xlabel('E[rDZ-rSib]');
    ylabel('frequency');
    hold on;
    yLimits = get(gca,'YLim');
    line([CorrDiff_DZ_Sib(N) CorrDiff_DZ_Sib(N)],[0 yLimits(2)],'Marker','.','Color','green');
    line([ctl_val_CorrDiff_DZ_Sib ctl_val_CorrDiff_DZ_Sib],[0 yLimits(2)],'Marker','.','LineStyle','-.','Color','red');
    hold off;
    set(gcf,'PaperPosition',[0 0 10 6])
    set(gcf,'PaperSize',[10 6])
    print('-dpdf',fullfile(ACEfit_Par.ResDir,['Diff_rDZ_rSib' ResSuf '.pdf']));
    
    
    Pvals_MZDZ  = [p_t_stats_MZ_DZ  p_CorrDiff_MZ_DZ]';
    Pvals_DZSib = [p_t_stats_DZ_Sib p_CorrDiff_DZ_Sib]';
    
    save(fullfile(ACEfit_Par.ResDir,['Pvals_AgHe' ResSuf]),'Pvals_MZDZ','Pvals_DZSib');
    
end


%
% Bootstrapping
%

if ACEfit_Par.nBoot>0
    
    nB = ACEfit_Par.nBoot+1;
    
    Boot_rMZ1  = zeros(nB,nMZF);
    Boot_rDZ1  = zeros(nB,nDZF);
    
    Boot_rDZ2  = zeros(nB,nDZF);
    Boot_rSib2 = zeros(nB,nSibP);
    
    Boot_rMZ1(end,:)  = rMZo;
    Boot_rDZ1(end,:)  = rDZo;
    
    Boot_rDZ2(end,:)  = rDZo;
    Boot_rSib2(end,:) = rSibo;
    
    % Compute rMZ & rDZ for all bootstrap replicates
    Boot_MZDZ = ACEfit_Par.Boot_index(:,1:nMZF+nDZF);
    for i=1:ACEfit_Par.nBoot
        Boot_label      = Boot_MZDZ(i,:);
        Boot_rMZ1(i,:)  = Boot_rMZ1(end,Boot_label(1:nMZF));
        Boot_rDZ1(i,:)  = Boot_rDZ1(end,Boot_label(nMZF+[1:nDZF])-nMZF);
    end
    
    % Compute rDZ & rSib for all bootstrap replicates
    for i=1:ACEfit_Par.nBoot
        Boot_label      = Boot_DZSib(i,:);
        Boot_rDZ2(i,:)  = Boot_rDZ2(end,Boot_label(1:nDZF));
        Boot_rSib2(i,:) = Boot_rSib2(end,Boot_label(nDZF+[1:nSibP])-nDZF);
    end
    
    % Bootstrapping mean differences
    Boot_diff_MZ_DZ  = mean(Boot_rMZ1,2)-mean(Boot_rDZ1,2);
    Boot_diff_DZ_Sib = mean(Boot_rDZ2,2)-mean(Boot_rSib2,2);
    
    % Confidence level 1-Balpha
    CI_MZDZ  = zeros(1,2);
    CI_DZSib = zeros(1,2);
    
    [CI_MZDZ(1),  CI_MZDZ(2) ] = Bootstrap_CI(Boot_diff_MZ_DZ,  Balpha);
    [CI_DZSib(1), CI_DZSib(2)] = Bootstrap_CI(Boot_diff_DZ_Sib, Balpha);
    
    save(fullfile(ACEfit_Par.ResDir,['CIs_AgHe' ResSuf]),'CI_MZDZ','CI_DZSib');
    
end


return



