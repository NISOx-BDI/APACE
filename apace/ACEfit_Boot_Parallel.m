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
    
    switch upper(ACEfit_Par.Model)
        
        case 'ACE'
            
            [MEANH2,WH2,MEDH2,Q3H2,MGMEDH2,MGQ3H2] = deal(zeros(nBootPerRun,1));
            [MEANC2,WC2,MEDC2,Q3C2,MGMEDC2,MGQ3C2] = deal(zeros(nBootPerRun,1));
            [MEANE2,WE2,MEDE2,Q3E2,MGMEDE2,MGQ3E2] = deal(zeros(nBootPerRun,1));
            
            for j=1:nBootPerRun
                
                if rem(j,50)==0; fprintf('%d ',j); end
                
                rand_pair = index(j,:);
                [Boot_label, Boot_kin] = deal([]);
                for nF=1:length(rand_pair)
                    IndF       = cFam(rand_pair(nF))+1:cFam(rand_pair(nF)+1);
                    Boot_label = [Boot_label IndF];
                    Boot_kin   = blkdiag(Boot_kin, ACEfit_Par.kin(IndF,IndF));
                end
                
                [SummaryACE] = ACEfit_Boot(ACEfit_Par,Boot_label,Boot_kin);
                
                MEANH2(j)  = SummaryACE(1,1);
                WH2(j)     = SummaryACE(2,1);
                MEDH2(j)   = SummaryACE(3,1);
                Q3H2(j)    = SummaryACE(4,1);
                MGMEDH2(j) = SummaryACE(5,1);
                MGQ3H2(j)  = SummaryACE(6,1);
                
                MEANC2(j)  = SummaryACE(1,2);
                WC2(j)     = SummaryACE(2,2);
                MEDC2(j)   = SummaryACE(3,2);
                Q3C2(j)    = SummaryACE(4,2);
                MGMEDC2(j) = SummaryACE(5,2);
                MGQ3C2(j)  = SummaryACE(6,2);
                
                MEANE2(j)  = SummaryACE(1,3);
                WE2(j)     = SummaryACE(2,3);
                MEDE2(j)   = SummaryACE(3,3);
                Q3E2(j)    = SummaryACE(4,3);
                MGMEDE2(j) = SummaryACE(5,3);
                MGQ3E2(j)  = SummaryACE(6,3);
                
            end
            
            str = fullfile(ACEfit_Par.ResDir,'BootCI_Parallel');
            save(sprintf('%s_%04d',str,k),'nBootPerRun','MEANH2','WH2','MEDH2','Q3H2','MGMEDH2','MGQ3H2',...
                                          'MEANC2','WC2','MEDC2','Q3C2','MGMEDC2','MGQ3C2',...
                                          'MEANE2','WE2','MEDE2','Q3E2','MGMEDE2','MGQ3E2');
            
        case 'AE'
            
            [MEANH2,WH2,MEDH2,Q3H2,MGMEDH2,MGQ3H2] = deal(zeros(nBootPerRun,1));
            [MEANE2,WE2,MEDE2,Q3E2,MGMEDE2,MGQ3E2] = deal(zeros(nBootPerRun,1));
            
            for j=1:nBootPerRun
                
                if rem(j,10)==0
                    fprintf('%d ',j);
                end
                % fprintf('%d',j);
                
                rand_pair = index(j,:);
                [Boot_label, Boot_kin] = deal([]);
                for nF=1:length(rand_pair)
                    IndF       = cFam(rand_pair(nF))+1:cFam(rand_pair(nF)+1);
                    Boot_label = [Boot_label IndF];
                    Boot_kin   = blkdiag(Boot_kin, ACEfit_Par.kin(IndF,IndF));
                end
                
                [SummaryAE] = ACEfit_Boot(ACEfit_Par,Boot_label,Boot_kin);
                
                MEANH2(j)  = SummaryAE(1,1);
                WH2(j)     = SummaryAE(2,1);
                MEDH2(j)   = SummaryAE(3,1);
                Q3H2(j)    = SummaryAE(4,1);
                MGMEDH2(j) = SummaryAE(5,1);
                MGQ3H2(j)  = SummaryAE(6,1);
                
                MEANE2(j)  = SummaryAE(1,2);
                WE2(j)     = SummaryAE(2,2);
                MEDE2(j)   = SummaryAE(3,2);
                Q3E2(j)    = SummaryAE(4,2);
                MGMEDE2(j) = SummaryAE(5,2);
                MGQ3E2(j)  = SummaryAE(6,2);
                
            end
            
            str = fullfile(ACEfit_Par.ResDir,'BootCI_Parallel');
            save(sprintf('%s_%04d',str,k),'nBootPerRun','MEANH2','WH2','MEDH2','Q3H2','MGMEDH2','MGQ3H2',...
                                                        'MEANE2','WE2','MEDE2','Q3E2','MGMEDE2','MGQ3E2');
            
            
    end
    
end

fprintf('\n');
return
