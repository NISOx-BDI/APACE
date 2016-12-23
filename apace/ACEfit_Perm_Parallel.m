function ACEfit_Perm_Parallel(ACEfit_Par,RunID)
%
% Parallel Computing.
%

nPm = ACEfit_Par.nPermPerRun;

fprintf('Permutation: ')

for k = RunID
    
    i0          = sum(nPm(1:(k-1)));
    nPermPerRun = nPm(k);
    index       = ACEfit_Par.Perm_index((i0+1):(i0+nPermPerRun),:);
    
    [MEANH2,WH2,MEDH2,Q3H2,MGMEDH2,MGQ3H2] = deal(zeros(nPermPerRun,1));
    
    if ~ACEfit_Par.NoImg
        uCnt = 0;
        maxT = zeros(nPermPerRun,1);
        if ACEfit_Par.Vs.ClustInf
            [maxK,maxM] = deal(zeros(nPermPerRun,1));
        end
    end
    
    for j = 1:nPermPerRun
        
        if rem(j,50)==0; fprintf('%d ',j); end
        
        Perm_label = index(j,:);
        
        [mK,mM,mT,uCount,SummaryA] = ACEfit_Perm(ACEfit_Par,Perm_label);
        
        if ~ACEfit_Par.NoImg
            maxT(j) = mT;
            uCnt    = uCnt + uCount;
            if ACEfit_Par.Vs.ClustInf
                maxK(j) = mK;
                maxM(j) = mM;
            end
        end
        
        MEANH2(j)  = SummaryA(1);
        WH2(j)     = SummaryA(2);
        MEDH2(j)   = SummaryA(3);
        Q3H2(j)    = SummaryA(4);
        MGMEDH2(j) = SummaryA(5);
        MGQ3H2(j)  = SummaryA(6);
        
    end
    
    str = fullfile(ACEfit_Par.ResDir,'ACEfit_Parallel');
    if ~ACEfit_Par.NoImg
        if ACEfit_Par.Vs.ClustInf
            save(sprintf('%s_%04d',str,k),'nPermPerRun','maxK','maxM','maxT','uCnt','MEANH2','WH2','MEDH2','Q3H2','MGMEDH2','MGQ3H2');
        else
            save(sprintf('%s_%04d',str,k),'nPermPerRun','maxT','uCnt','MEANH2','WH2','MEDH2','Q3H2','MGMEDH2','MGQ3H2');
        end
    else
        save(sprintf('%s_%04d',str,k),'nPermPerRun','MEANH2','WH2','MEDH2','Q3H2','MGMEDH2','MGQ3H2');
    end
    
end

fprintf('\n');
return
