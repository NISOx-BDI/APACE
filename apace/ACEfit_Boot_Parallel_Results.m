function ACEfit_Boot_Parallel_Results(ACEfit_Par)
%
% Merge together results of RunID.
%

switch upper(ACEfit_Par.Model)
    
    case 'ACE'
        
        mean_A  = zeros(ACEfit_Par.nBoot,1);
        wh2_A   = zeros(ACEfit_Par.nBoot,1);
        med_A   = zeros(ACEfit_Par.nBoot,1);
        q3_A    = zeros(ACEfit_Par.nBoot,1);
        mGmed_A = zeros(ACEfit_Par.nBoot,1);
        mGq3_A  = zeros(ACEfit_Par.nBoot,1);
        
        mean_C  = zeros(ACEfit_Par.nBoot,1);
        wc2_C   = zeros(ACEfit_Par.nBoot,1);
        med_C   = zeros(ACEfit_Par.nBoot,1);
        q3_C    = zeros(ACEfit_Par.nBoot,1);
        mGmed_C = zeros(ACEfit_Par.nBoot,1);
        mGq3_C  = zeros(ACEfit_Par.nBoot,1);
        
        mean_E  = zeros(ACEfit_Par.nBoot,1);
        we2_E   = zeros(ACEfit_Par.nBoot,1);
        med_E   = zeros(ACEfit_Par.nBoot,1);
        q3_E    = zeros(ACEfit_Par.nBoot,1);
        mGmed_E = zeros(ACEfit_Par.nBoot,1);
        mGq3_E  = zeros(ACEfit_Par.nBoot,1);
        
        i0 = 0;
        for i=ACEfit_Par.RunID
            try
                file = fullfile(ACEfit_Par.ResDir,sprintf('BootCI_Parallel_%04d',i));
            catch
                file = fullfile(ACEfit_Par.ResDir,sprintf('BootCI_Parallel_%4d',i));
            end
            s = load(file);
            
            mean_A( (i0+1):(i0+s.nBootPerRun)) = s.MEANH2;
            wh2_A(  (i0+1):(i0+s.nBootPerRun)) = s.WH2;
            med_A(  (i0+1):(i0+s.nBootPerRun)) = s.MEDH2;
            q3_A(   (i0+1):(i0+s.nBootPerRun)) = s.Q3H2;
            mGmed_A((i0+1):(i0+s.nBootPerRun)) = s.MGMEDH2;
            mGq3_A( (i0+1):(i0+s.nBootPerRun)) = s.MGQ3H2;
            
            mean_C( (i0+1):(i0+s.nBootPerRun)) = s.MEANC2;
            wc2_C(  (i0+1):(i0+s.nBootPerRun)) = s.WC2;
            med_C(  (i0+1):(i0+s.nBootPerRun)) = s.MEDC2;
            q3_C(   (i0+1):(i0+s.nBootPerRun)) = s.Q3C2;
            mGmed_C((i0+1):(i0+s.nBootPerRun)) = s.MGMEDC2;
            mGq3_C( (i0+1):(i0+s.nBootPerRun)) = s.MGQ3C2;
            
            mean_E( (i0+1):(i0+s.nBootPerRun)) = s.MEANE2;
            we2_E(  (i0+1):(i0+s.nBootPerRun)) = s.WE2;
            med_E(  (i0+1):(i0+s.nBootPerRun)) = s.MEDE2;
            q3_E(   (i0+1):(i0+s.nBootPerRun)) = s.Q3E2;
            mGmed_E((i0+1):(i0+s.nBootPerRun)) = s.MGMEDE2;
            mGq3_E( (i0+1):(i0+s.nBootPerRun)) = s.MGQ3E2;
            
            i0 = i0 + s.nBootPerRun;
        end
        clear s
        
        meanh2_ACE  = [mean_A;  ACEfit_Par.SummaryA(1)];
        wh2_ACE     = [wh2_A;   ACEfit_Par.SummaryA(2)];
        medh2_ACE   = [med_A;   ACEfit_Par.SummaryA(3)];
        q3h2_ACE    = [q3_A;    ACEfit_Par.SummaryA(4)];
        mGmedh2_ACE = [mGmed_A; ACEfit_Par.SummaryA(5)];
        mGq3h2_ACE  = [mGq3_A;  ACEfit_Par.SummaryA(6)];
        
        meanc2_ACE  = [mean_C;  ACEfit_Par.SummaryC(1)];
        wc2_ACE     = [wc2_C;   ACEfit_Par.SummaryC(2)];
        medc2_ACE   = [med_C;   ACEfit_Par.SummaryC(3)];
        q3c2_ACE    = [q3_C;    ACEfit_Par.SummaryC(4)];
        mGmedc2_ACE = [mGmed_C; ACEfit_Par.SummaryC(5)];
        mGq3c2_ACE  = [mGq3_C;  ACEfit_Par.SummaryC(6)];
        
        meane2_ACE  = [mean_E;  ACEfit_Par.SummaryE(1)];
        we2_ACE     = [we2_E;   ACEfit_Par.SummaryE(2)];
        mede2_ACE   = [med_E;   ACEfit_Par.SummaryE(3)];
        q3e2_ACE    = [q3_E;    ACEfit_Par.SummaryE(4)];
        mGmede2_ACE = [mGmed_E; ACEfit_Par.SummaryE(5)];
        mGq3e2_ACE  = [mGq3_E;  ACEfit_Par.SummaryE(6)];
        
        save(fullfile(ACEfit_Par.ResDir,'ACEfit_Boot'),'meanh2_ACE','wh2_ACE','medh2_ACE','q3h2_ACE','mGmedh2_ACE','mGq3h2_ACE',...
                                                       'meanc2_ACE','wc2_ACE','medc2_ACE','q3c2_ACE','mGmedc2_ACE','mGq3c2_ACE',...
                                                       'meane2_ACE','we2_ACE','mede2_ACE','q3e2_ACE','mGmede2_ACE','mGq3e2_ACE');
        
    case 'AE'
        
        mean_A  = zeros(ACEfit_Par.nBoot,1);
        wh2_A   = zeros(ACEfit_Par.nBoot,1);
        med_A   = zeros(ACEfit_Par.nBoot,1);
        q3_A    = zeros(ACEfit_Par.nBoot,1);
        mGmed_A = zeros(ACEfit_Par.nBoot,1);
        mGq3_A  = zeros(ACEfit_Par.nBoot,1);
        
        mean_E  = zeros(ACEfit_Par.nBoot,1);
        we2_E   = zeros(ACEfit_Par.nBoot,1);
        med_E   = zeros(ACEfit_Par.nBoot,1);
        q3_E    = zeros(ACEfit_Par.nBoot,1);
        mGmed_E = zeros(ACEfit_Par.nBoot,1);
        mGq3_E  = zeros(ACEfit_Par.nBoot,1);
        
        i0 = 0;
        for i=ACEfit_Par.RunID
            try
                file = fullfile(ACEfit_Par.ResDir,sprintf('BootCI_Parallel_%04d',i));
            catch
                file = fullfile(ACEfit_Par.ResDir,sprintf('BootCI_Parallel_%4d',i));
            end
            s = load(file);
            
            mean_A( (i0+1):(i0+s.nBootPerRun)) = s.MEANH2;
            wh2_A(  (i0+1):(i0+s.nBootPerRun)) = s.WH2;
            med_A(  (i0+1):(i0+s.nBootPerRun)) = s.MEDH2;
            q3_A(   (i0+1):(i0+s.nBootPerRun)) = s.Q3H2;
            mGmed_A((i0+1):(i0+s.nBootPerRun)) = s.MGMEDH2;
            mGq3_A( (i0+1):(i0+s.nBootPerRun)) = s.MGQ3H2;
            
            mean_E( (i0+1):(i0+s.nBootPerRun)) = s.MEANE2;
            we2_E(  (i0+1):(i0+s.nBootPerRun)) = s.WE2;
            med_E(  (i0+1):(i0+s.nBootPerRun)) = s.MEDE2;
            q3_E(   (i0+1):(i0+s.nBootPerRun)) = s.Q3E2;
            mGmed_E((i0+1):(i0+s.nBootPerRun)) = s.MGMEDE2;
            mGq3_E( (i0+1):(i0+s.nBootPerRun)) = s.MGQ3E2;
            
            i0 = i0 + s.nBootPerRun;
        end
        clear s
        
        meanh2_ACE  = [mean_A;  ACEfit_Par.SummaryA(1)];
        wh2_ACE     = [wh2_A;   ACEfit_Par.SummaryA(2)];
        medh2_ACE   = [med_A;   ACEfit_Par.SummaryA(3)];
        q3h2_ACE    = [q3_A;    ACEfit_Par.SummaryA(4)];
        mGmedh2_ACE = [mGmed_A; ACEfit_Par.SummaryA(5)];
        mGq3h2_ACE  = [mGq3_A;  ACEfit_Par.SummaryA(6)];
        
        meane2_ACE  = [mean_E;  ACEfit_Par.SummaryE(1)];
        we2_ACE     = [we2_E;   ACEfit_Par.SummaryE(2)];
        mede2_ACE   = [med_E;   ACEfit_Par.SummaryE(3)];
        q3e2_ACE    = [q3_E;    ACEfit_Par.SummaryE(4)];
        mGmede2_ACE = [mGmed_E; ACEfit_Par.SummaryE(5)];
        mGq3e2_ACE  = [mGq3_E;  ACEfit_Par.SummaryE(6)];
        
        save(fullfile(ACEfit_Par.ResDir,'ACEfit_Boot'),'meanh2_ACE','wh2_ACE','medh2_ACE','q3h2_ACE','mGmedh2_ACE','mGq3h2_ACE',...
                                                       'meane2_ACE','we2_ACE','mede2_ACE','q3e2_ACE','mGmede2_ACE','mGq3e2_ACE');
        
end

return
