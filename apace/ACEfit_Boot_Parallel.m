function ACEfit_Boot_Parallel(ACEfit_Par,RunID)
%
% Parallel Computation.
%
%_______________________________________________________________________
% Version: http://github.com/NISOx-BDI/APACE/tree/$Format:%h$
%          $Format:%ci$

nBt = ACEfit_Par.nBootPerRun;

IndFam_end = cumsum(ACEfit_Par.nFam);
IndFam_beg = IndFam_end-ACEfit_Par.nFam+1;

fprintf('Bootstrap: ')

for k = RunID
    
    i0          = sum(nBt(1:(k-1)));
    nBootPerRun = nBt(k);
    index       = ACEfit_Par.Boot_index((i0+1):(i0+nBootPerRun),:);
    
    [MEANH2,WH2,MEDH2,Q3H2,MGMEDH2,MGQ3H2] = deal(zeros(nBootPerRun,1));
    [MEANC2,WC2,MEDC2,Q3C2,MGMEDC2,MGQ3C2] = deal(zeros(nBootPerRun,1));
    [MEANE2,WE2,MEDE2,Q3E2,MGMEDE2,MGQ3E2] = deal(zeros(nBootPerRun,1));
    
    for j = 1:nBootPerRun
        
        if ~rem(j,50); fprintf('%d ',j); end
        
        % Create kinship matrix for each bootstrap replicate
        rand_pair = index(j,:);
        
        BootFam_end = IndFam_end(rand_pair);
        BootFam_beg = IndFam_beg(rand_pair);
        
        Boot_nFam = ACEfit_Par.nFam(rand_pair);
        nB        = sum(Boot_nFam);
        
        Boot_label = zeros(1,nB);
        Boot_kin   = zeros(nB);
        iSubj      = 0;
        for nF = 1:length(rand_pair)
           IndexF = BootFam_beg(nF):BootFam_end(nF); 
           
           Boot_label(iSubj+[1:Boot_nFam(nF)])                         = IndexF;
           Boot_kin(  iSubj+[1:Boot_nFam(nF)],iSubj+[1:Boot_nFam(nF)]) = ACEfit_Par.kin(IndexF,IndexF);
           
           iSubj = iSubj + Boot_nFam(nF);
        end

        [SummaryACE] = ACEfit_Boot(ACEfit_Par,Boot_label,Boot_kin);
        
        switch upper(ACEfit_Par.Model)
            
            case 'ACE'
                
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
                
            case 'AE'
                
                MEANH2(j)  = SummaryACE(1,1);
                WH2(j)     = SummaryACE(2,1);
                MEDH2(j)   = SummaryACE(3,1);
                Q3H2(j)    = SummaryACE(4,1);
                MGMEDH2(j) = SummaryACE(5,1);
                MGQ3H2(j)  = SummaryACE(6,1);
                
                MEANE2(j)  = SummaryACE(1,2);
                WE2(j)     = SummaryACE(2,2);
                MEDE2(j)   = SummaryACE(3,2);
                Q3E2(j)    = SummaryACE(4,2);
                MGMEDE2(j) = SummaryACE(5,2);
                MGQ3E2(j)  = SummaryACE(6,2);
                
        end
        
    end
    
    
    switch upper(ACEfit_Par.Model)
        case 'ACE'
            str = fullfile(ACEfit_Par.ResDir,'BootCI_Parallel');
            save(sprintf('%s_%04d',str,k),'nBootPerRun','MEANH2','WH2','MEDH2','Q3H2','MGMEDH2','MGQ3H2',...
                                                        'MEANC2','WC2','MEDC2','Q3C2','MGMEDC2','MGQ3C2',...
                                                        'MEANE2','WE2','MEDE2','Q3E2','MGMEDE2','MGQ3E2');
        case 'AE'
            str = fullfile(ACEfit_Par.ResDir,'BootCI_Parallel');
            save(sprintf('%s_%04d',str,k),'nBootPerRun','MEANH2','WH2','MEDH2','Q3H2','MGMEDH2','MGQ3H2',...
                                                        'MEANE2','WE2','MEDE2','Q3E2','MGMEDE2','MGQ3E2');
    end
    
end

fprintf('\n')
return
