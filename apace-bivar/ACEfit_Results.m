function ACEfit_Results(ACEfit_Par,FWEalpha,FDRalpha)
%
% Generate Permutation Results
%
I_data = ACEfit_Par.I_data;

load(fullfile(ACEfit_Par.ResDir,'ACEfit_Perm'));
N             = ACEfit_Par.nPerm+1;
ctl_val_index = N-floor(N*FWEalpha);

f = 0;

%
% (1) Maximum Test Statistic
%
f = f+1;
figure(f);
[sT,oT]   = sort(max_T_ERV);
ctl_val_T = sT(ctl_val_index);
[f1,x1]   = hist(sT,50);
n_T       = min(find(sT(:)==max_T_ERV(N)));
bar(x1,f1/N);

% Calculate the permutation-based p-value
p_T = (N-n_T+1)/N;

set(gca, 'xtick', 0)
if ctl_val_T==max_T_ERV(N)
    set(gca, 'xtick', ctl_val_T)
else
    set(gca, 'xtick', sort([ctl_val_T max_T_ERV(N)]))
end

title(sprintf('Empirical Distribution of Maximum Test Statistic of ERV, FWE P-value = %.3f',p_T));
xlabel('T');
ylabel('f(T)');
hold on;
yLimits = get(gca,'YLim');
% Plot the original point (green)
line([max_T_ERV(N) max_T_ERV(N)],[0 yLimits(2)],'Marker','.','Color','green');
% Plot the critical threshold (red)
line([ctl_val_T ctl_val_T],[0 yLimits(2)],'Marker','.','LineStyle','-.','Color','red');
hold off;
set(gcf,'PaperPosition',[0 0 10 6])
set(gcf,'PaperSize',[10 6])
print('-dpdf',fullfile(ACEfit_Par.ResDir,'maxT_ERV.pdf'));

% Voxel-wise FDR correction
Pval   = unPval_ERV(I_data,I_data);
Pid    = find(triu(ones(length(I_data)),1));
Pval   = Pval(Pid);
unPval = Pval;
nMsk   = numel(Pval);

% Compute FDR corrected p-values
[sP,iP]  = sort(Pval);
cV       = 1;
Qs       = min(sP./(1:nMsk)'*nMsk*cV,1);
sP(end)  = Qs(end);
for iMsk = nMsk-1:-1:1
    sP(iMsk) = min(Qs(iMsk),sP(iMsk+1));
end

FDR_min  = min(Qs);
Pval(iP) = sP;

%%% FDR plot
f = f+1;
figure(f);
plot((1:nMsk)'/nMsk,sort(unPval),(1:nMsk)'/nMsk,sort(Pval));
legend({'Uncorr P','FDR-corr P'},'Location','NorthWest')
xlabel('Expected Uncorrected Ordered P');
ylabel('Ordered P');
title('Element-wise P-values for ERV: FDR plot');
axis equal;axis tight;axis([0 1 0 1]);
abline(0,1,'LineStyle','-','color',[.5 .5 .5]);
abline(0,FDRalpha,'LineStyle',':','color','red')
set(gcf,'PaperPosition',[0 0 6 6])
set(gcf,'PaperSize',[6 6])
print('-dpdf',fullfile(ACEfit_Par.ResDir,'FDRplot_ERV.pdf'));

FDR_P      = ones(length(I_data));
FDR_P(Pid) = Pval;

FDR_Pval   = ones(size(unPval_ERV));
FDR_Pval(I_data,I_data) = FDR_P;

% -log10(FDR corrected p-value)
FDR_Pval_ERV = -log10(FDR_Pval);

fprintf('The best FDR attainable is %.2f. \n', FDR_min);

% -log10(uncorrected P-value)
unPval_ERV   = -log10(unPval_ERV);

% Voxel-wise FWE correction
Tstat        = ACEfit_Par.Stats;
corrPval_ERV = zeros(size(unPval_ERV));

for id = 1:length(I_data)
    i = I_data(id);
    for j = I_data(id+1:end)
        % -log10(FWE corrected P-value)
        corrPval_ERV(i,j) = -log10(mean(Tstat(i,j)<=max_T_ERV));
    end
end


%
% Save and output the permutation results
%
save(fullfile(ACEfit_Par.ResDir,'Perm_Results'),'unPval_ERV','corrPval_ERV','FDR_Pval_ERV','-v7.3');


return
