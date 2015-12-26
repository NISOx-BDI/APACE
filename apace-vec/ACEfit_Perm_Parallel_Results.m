function ACEfit_Perm_Parallel_Results(ACEfit_Par)
%
% Merge together results of RunID.
%

mean_A  = zeros(ACEfit_Par.nPerm,1);
wh2_A   = zeros(ACEfit_Par.nPerm,1);
med_A   = zeros(ACEfit_Par.nPerm,1);
q3_A    = zeros(ACEfit_Par.nPerm,1);
mGmed_A = zeros(ACEfit_Par.nPerm,1);
mGq3_A  = zeros(ACEfit_Par.nPerm,1);

if ~ACEfit_Par.NoImg
    max_T = zeros(ACEfit_Par.nPerm,1);
    Count = 0;
    if ACEfit_Par.Vs.ClustInf
        max_K = zeros(ACEfit_Par.nPerm,1);
        max_M = zeros(ACEfit_Par.nPerm,1);
    end
end

i0 = 0;
for i=ACEfit_Par.RunID
    try
        file = fullfile(ACEfit_Par.ResDir,sprintf('ACEfit_Parallel_%04d',i));
    catch
        file = fullfile(ACEfit_Par.ResDir,sprintf('ACEfit_Parallel_%4d',i));
    end
    s    = load(file);
    
    mean_A( (i0+1):(i0+s.nPermPerRun)) = s.MEANH2;
    wh2_A(  (i0+1):(i0+s.nPermPerRun)) = s.WH2;
    med_A(  (i0+1):(i0+s.nPermPerRun)) = s.MEDH2;
    q3_A(   (i0+1):(i0+s.nPermPerRun)) = s.Q3H2;
    mGmed_A((i0+1):(i0+s.nPermPerRun)) = s.MGMEDH2;
    mGq3_A( (i0+1):(i0+s.nPermPerRun)) = s.MGQ3H2;
    
    if ~ACEfit_Par.NoImg
        max_T((i0+1):(i0+s.nPermPerRun)) = s.maxT;
        Count                            = Count + s.uCnt;
        if ACEfit_Par.Vs.ClustInf
            max_K((i0+1):(i0+s.nPermPerRun)) = s.maxK;
            max_M((i0+1):(i0+s.nPermPerRun)) = s.maxM;
        end
    end
    
    i0   = i0 + s.nPermPerRun;
end
clear s

mean_ACE  = [mean_A;  ACEfit_Par.SummaryA(1)];
wh2_ACE   = [wh2_A;   ACEfit_Par.SummaryA(2)];
med_ACE   = [med_A;   ACEfit_Par.SummaryA(3)];
q3_ACE    = [q3_A;    ACEfit_Par.SummaryA(4)];
mGmed_ACE = [mGmed_A; ACEfit_Par.SummaryA(5)];
mGq3_ACE  = [mGq3_A;  ACEfit_Par.SummaryA(6)];

if ~ACEfit_Par.NoImg
    
    max_T_ACE  = [max_T; ACEfit_Par.mT];
    
    % Compute the uncorrected empirical p-value
    unPval_ACE = (Count+1)/(ACEfit_Par.nPerm+1);
    
    if ACEfit_Par.Vs.ClustInf

        max_K_ACE = [max_K; ACEfit_Par.mK];
        max_M_ACE = [max_M; ACEfit_Par.mM];
        
        save(fullfile(ACEfit_Par.ResDir,'ACEfit_Perm'),'max_K_ACE','max_M_ACE','max_T_ACE','unPval_ACE',...
                                                       'mean_ACE','wh2_ACE','med_ACE','q3_ACE','mGmed_ACE','mGq3_ACE');
    else
        save(fullfile(ACEfit_Par.ResDir,'ACEfit_Perm'),'max_T_ACE','unPval_ACE',...
                                                       'mean_ACE','wh2_ACE','med_ACE','q3_ACE','mGmed_ACE','mGq3_ACE');
    end
    
else
    save(fullfile(ACEfit_Par.ResDir,'ACEfit_Perm'),'mean_ACE','wh2_ACE','med_ACE','q3_ACE','mGmed_ACE','mGq3_ACE');
    
end

return
