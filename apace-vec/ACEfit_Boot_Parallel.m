function ACEfit_Boot_Parallel(ACEfit_Par,RunID)
%
% Parallel Computation.
%

cFam = [0; cumsum(ACEfit_Par.nFam)];
nBt  = ACEfit_Par.nBootPerRun;

for k=RunID
    
    i0          = sum(nBt(1:(k-1)));
    nBootPerRun = nBt(k);
    index       = ACEfit_Par.Boot_index((i0+1):(i0+nBootPerRun),:);
    
    [MEANH2,WH2,MEDH2,Q3H2,MGMEDH2,MGQ3H2] = deal(zeros(nBootPerRun,1));
    [MEANC2,WC2,MEDC2,Q3C2,MGMEDC2,MGQ3C2] = deal(zeros(nBootPerRun,1));
    [MEANE2,WE2,MEDE2,Q3E2,MGMEDE2,MGQ3E2] = deal(zeros(nBootPerRun,1));
    
    for j=1:nBootPerRun
        
        if ~rem(j,50); fprintf('%d ',j); end
        
        rand_pair = index(j,:);
        [Boot_label, Boot_kin] = deal([]);
        for nF=1:length(rand_pair)
            IndF       = cFam(rand_pair(nF))+1:cFam(rand_pair(nF)+1);
            Boot_label = [Boot_label IndF];
            Boot_kin   = blkdiag(Boot_kin, ACEfit_Par.kin(IndF,IndF));
        end
        
        [SummaryA,SummaryC,SummaryE] = ACEfit_Boot(ACEfit_Par,Boot_label,Boot_kin);
        
        MEANH2(j)  = SummaryA(1);
        WH2(j)     = SummaryA(2);
        MEDH2(j)   = SummaryA(3);
        Q3H2(j)    = SummaryA(4);
        MGMEDH2(j) = SummaryA(5);
        MGQ3H2(j)  = SummaryA(6);
        
        MEANC2(j)  = SummaryC(1);
        WC2(j)     = SummaryC(2);
        MEDC2(j)   = SummaryC(3);
        Q3C2(j)    = SummaryC(4);
        MGMEDC2(j) = SummaryC(5);
        MGQ3C2(j)  = SummaryC(6);
        
        MEANE2(j)  = SummaryE(1);
        WE2(j)     = SummaryE(2);
        MEDE2(j)   = SummaryE(3);
        Q3E2(j)    = SummaryE(4);
        MGMEDE2(j) = SummaryE(5);
        MGQ3E2(j)  = SummaryE(6);
        
    end
    
    fprintf('\n')
    
    str = fullfile(ACEfit_Par.ResDir,'BootCI_Parallel');
    save(sprintf('%s_%04d',str,k),'nBootPerRun','MEANH2','WH2','MEDH2','Q3H2','MGMEDH2','MGQ3H2',...
                                                'MEANC2','WC2','MEDC2','Q3C2','MGMEDC2','MGQ3C2',...
                                                'MEANE2','WE2','MEDE2','Q3E2','MGMEDE2','MGQ3E2');
    
end

return
