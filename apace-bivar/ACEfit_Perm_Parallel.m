function ACEfit_Perm_Parallel(ACEfit_Par,RunID)
%
% Parallel Computing.
%

nPm = ACEfit_Par.nPermPerRun;

for k = RunID
    
    nPermPerRun = nPm(k);
    
    i0    = sum(nPm(1:k-1));
    index = ACEfit_Par.Perm_index(i0+1:i0+nPermPerRun,:);
    
    uCnt = 0;
    maxT = zeros(nPermPerRun,1);
    
    for j = 1:nPermPerRun
        
        fprintf('%d',j);
        Perm_label = index(j,:);
        
        [mT,uCount] = ACEfit_Perm(ACEfit_Par,Perm_label);
        
        maxT(j) = mT;
        uCnt    = uCnt + uCount;
        
    end
    
    str = fullfile(ACEfit_Par.ResDir,'ACEfit_Parallel');    
    save(sprintf('%s_%04d',str,k),'nPermPerRun','maxT','uCnt','-v7.3');
    
end

return