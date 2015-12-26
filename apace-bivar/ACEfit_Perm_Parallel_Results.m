function ACEfit_Perm_Parallel_Results(ACEfit_Par)
%
% Merge together results of RunID.
%

Count = 0;
max_T = zeros(ACEfit_Par.nPerm,1);

i0 = 0;
for i = ACEfit_Par.RunID
    
    s = load(fullfile(ACEfit_Par.ResDir,sprintf('ACEfit_Parallel_%04d',i)));
    
    max_T(i0+[1:s.nPermPerRun]) = s.maxT;
    Count                       = Count + s.uCnt;
    
    i0 = i0 + s.nPermPerRun;
    
end
clear s

% Compute the maximum test statistic & uncorrected empirical p-value
max_T_ERV  = [max_T; ACEfit_Par.mT];
unPval_ERV = (Count+1)/(ACEfit_Par.nPerm+1);

save(fullfile(ACEfit_Par.ResDir,'ACEfit_Perm'),'max_T_ERV','unPval_ERV','-v7.3');

return
