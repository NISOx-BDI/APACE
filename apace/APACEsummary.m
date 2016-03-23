function [varargout] = APACEsummary(ACEfit_Par,varargin)
%
% APACEsummary(ACEfit_Par,'ResultsFile')   - Saves results in ResultsFile.csv
% APACEsummary(ACEfit_Par)                 - Writes text results to stdout
%
% [Ests,CIs,Ps,Names] = APACEsummary(ACEfit_Par) - Also return values
%

if nargin<=1
    StdOut  = 1;
else
    StdOut  = 0;
    
    if strcmpi(spm_str_manip(varargin{1},'e'),'csv');
        CSVname = varargin{1};
    else
        CSVname = [varargin{1} '.csv'];
    end
    
    OutFile = fullfile(ACEfit_Par.ResDir,CSVname);
end

Names = {'Mean','VarWtMean','AgHe','Median Q2','Q3','Mean(>Q2)','Mean(>Q3)'};
% Names = {'Mean','VarWtMean','AgHeNorm','AgHeNoNorm','Median Q2','Q3','Mean(>Q2)','Mean(>Q3)'};

if exist(fullfile(ACEfit_Par.ResDir,'Ests_AgHe.mat'),'file')
    E_s = load(fullfile(ACEfit_Par.ResDir,'Ests_AgHe'));        % Est_MZDZ Est_DZSib
    % Ne = load(fullfile(ACEfit_Par.ResDir,'Ests_AgHe_Norm'));
    % Oe = load(fullfile(ACEfit_Par.ResDir,'Ests_AgHe_NoNorm'));
else
    E_s.Est_MZDZ  = NaN;
    E_s.Est_DZSib = NaN;
end

switch upper(ACEfit_Par.Model)
    
    case 'ACE'
        
        Ests = [ ACEfit_Par.SummaryA(1), ACEfit_Par.SummaryC(1), ACEfit_Par.SummaryE(1) ;
                 ACEfit_Par.SummaryA(2), ACEfit_Par.SummaryC(2), ACEfit_Par.SummaryE(2) ;
                 E_s.Est_MZDZ*2,         E_s.Est_DZSib,          NaN                    ;
                 ACEfit_Par.SummaryA(3), ACEfit_Par.SummaryC(3), ACEfit_Par.SummaryE(3) ;
                 ACEfit_Par.SummaryA(4), ACEfit_Par.SummaryC(4), ACEfit_Par.SummaryE(4) ;
                 ACEfit_Par.SummaryA(5), ACEfit_Par.SummaryC(5), ACEfit_Par.SummaryE(5) ;
                 ACEfit_Par.SummaryA(6), ACEfit_Par.SummaryC(6), ACEfit_Par.SummaryE(6) ];
        % Ests = [ ACEfit_Par.SummaryA(1), ACEfit_Par.SummaryC(1), ACEfit_Par.SummaryE(1) ;
        %          ACEfit_Par.SummaryA(2), ACEfit_Par.SummaryC(2), ACEfit_Par.SummaryE(2) ;
        %          Ne.Est_MZDZ*2,          Ne.Est_DZSib,           NaN                    ;
        %          Oe.Est_MZDZ*2,          Oe.Est_DZSib,           NaN                    ;
        %          ACEfit_Par.SummaryA(3), ACEfit_Par.SummaryC(3), ACEfit_Par.SummaryE(3) ;
        %          ACEfit_Par.SummaryA(4), ACEfit_Par.SummaryC(4), ACEfit_Par.SummaryE(4) ;
        %          ACEfit_Par.SummaryA(5), ACEfit_Par.SummaryC(5), ACEfit_Par.SummaryE(5) ;
        %          ACEfit_Par.SummaryA(6), ACEfit_Par.SummaryC(6), ACEfit_Par.SummaryE(6) ];
        
    case 'AE'
        
        Ests = [ ACEfit_Par.SummaryA(1), NaN,           ACEfit_Par.SummaryE(1) ;
                 ACEfit_Par.SummaryA(2), NaN,           ACEfit_Par.SummaryE(2) ;
                 E_s.Est_MZDZ*2,         E_s.Est_DZSib, NaN                    ;
                 ACEfit_Par.SummaryA(3), NaN,           ACEfit_Par.SummaryE(3) ;
                 ACEfit_Par.SummaryA(4), NaN,           ACEfit_Par.SummaryE(4) ;
                 ACEfit_Par.SummaryA(5), NaN,           ACEfit_Par.SummaryE(5) ;
                 ACEfit_Par.SummaryA(6), NaN,           ACEfit_Par.SummaryE(6) ];
        % Ests = [ ACEfit_Par.SummaryA(1), NaN,          ACEfit_Par.SummaryE(1) ;
        %          ACEfit_Par.SummaryA(2), NaN,          ACEfit_Par.SummaryE(2) ;
        %          Ne.Est_MZDZ*2,          Ne.Est_DZSib, NaN                    ;
        %          Oe.Est_MZDZ*2,          Oe.Est_DZSib, NaN                    ;
        %          ACEfit_Par.SummaryA(3), NaN,          ACEfit_Par.SummaryE(3) ;
        %          ACEfit_Par.SummaryA(4), NaN,          ACEfit_Par.SummaryE(4) ;
        %          ACEfit_Par.SummaryA(5), NaN,          ACEfit_Par.SummaryE(5) ;
        %          ACEfit_Par.SummaryA(6), NaN,          ACEfit_Par.SummaryE(6) ];
        
end

if ACEfit_Par.nPerm>0
    
    load(fullfile(ACEfit_Par.ResDir,'Pvals_h2'))               % Pvals_h2
    
    if exist(fullfile(ACEfit_Par.ResDir,'Pvals_AgHe.mat'),'file')
        P_s = load(fullfile(ACEfit_Par.ResDir,'Pvals_AgHe'));  % Pvals_MZDZ Pvals_DZSib
        % Np = load(fullfile(ACEfit_Par.ResDir,'Pvals_AgHe_Norm'));
        % Op = load(fullfile(ACEfit_Par.ResDir,'Pvals_AgHe_NoNorm'));
    else
        P_s.Pvals_MZDZ  = NaN;
        P_s.Pvals_DZSib = NaN;
    end  
    
    Ps = [ Pvals_h2(1),       NaN,                NaN ;
           Pvals_h2(2),       NaN,                NaN ;
           P_s.Pvals_MZDZ(1), P_s.Pvals_DZSib(1), NaN ;
           Pvals_h2(3),       NaN,                NaN ;
           Pvals_h2(4),       NaN,                NaN ;
           Pvals_h2(5),       NaN,                NaN ;
           Pvals_h2(6),       NaN,                NaN ];
    % Ps = [ Pvals_h2(1),      NaN,               NaN ;
    %        Pvals_h2(2),      NaN,               NaN ;
    %        Np.Pvals_MZDZ(1), Np.Pvals_DZSib(1), NaN ;
    %        Op.Pvals_MZDZ(1), Op.Pvals_DZSib(1), NaN ;
    %        Pvals_h2(3),      NaN,               NaN ;
    %        Pvals_h2(4),      NaN,               NaN ;
    %        Pvals_h2(5),      NaN,               NaN ;
    %        Pvals_h2(6),      NaN,               NaN ];
    
else
    
    Ps = NaN(length(Names),3);

end

if ACEfit_Par.nBoot>0
    
    load(fullfile(ACEfit_Par.ResDir,'Boot_CIs'))             % alpha CIs_h2 CIs_c2 CIs_e2
    
    if exist(fullfile(ACEfit_Par.ResDir,'CIs_AgHe.mat'),'file')
        CI_s = load(fullfile(ACEfit_Par.ResDir,'CIs_AgHe')); % CI_MZDZ CI_DZSib
        % Nc = load(fullfile(ACEfit_Par.ResDir,'CIs_AgHe_Norm'));
        % Oc = load(fullfile(ACEfit_Par.ResDir,'CIs_AgHe_NoNorm'));
    else
        CI_s.CI_MZDZ =  [NaN NaN];
        CI_s.CI_DZSib = [NaN NaN];
    end  
    
    switch upper(ACEfit_Par.Model)
        
        case 'ACE'
            
            CIs = { CIs_h2(1,:),    CIs_c2(1,:),   CIs_e2(1,:) ;
                    CIs_h2(2,:),    CIs_c2(2,:),   CIs_e2(2,:) ;
                    CI_s.CI_MZDZ*2, CI_s.CI_DZSib, [NaN NaN]   ;
                    CIs_h2(3,:),    CIs_c2(3,:),   CIs_e2(3,:) ;
                    CIs_h2(4,:),    CIs_c2(4,:),   CIs_e2(4,:) ;
                    CIs_h2(5,:),    CIs_c2(5,:),   CIs_e2(5,:) ;
                    CIs_h2(6,:),    CIs_c2(6,:),   CIs_e2(6,:) };
            % CIs = { CIs_h2(1,:),  CIs_c2(1,:), CIs_e2(1,:) ;
            %         CIs_h2(2,:),  CIs_c2(2,:), CIs_e2(2,:) ;
            %         Nc.CI_MZDZ*2, Nc.CI_DZSib, [NaN NaN]   ;
            %         Oc.CI_MZDZ*2, Oc.CI_DZSib, [NaN NaN]   ;
            %         CIs_h2(3,:),  CIs_c2(3,:), CIs_e2(3,:) ;
            %         CIs_h2(4,:),  CIs_c2(4,:), CIs_e2(4,:) ;
            %         CIs_h2(5,:),  CIs_c2(5,:), CIs_e2(5,:) ;
            %         CIs_h2(6,:),  CIs_c2(6,:), CIs_e2(6,:) };
            
        case 'AE'
            
            CIs = { CIs_h2(1,:),    [NaN NaN],     CIs_e2(1,:) ;
                    CIs_h2(2,:),    [NaN NaN],     CIs_e2(2,:) ;
                    CI_s.CI_MZDZ*2, CI_s.CI_DZSib, [NaN NaN]   ;
                    CIs_h2(3,:),    [NaN NaN],     CIs_e2(3,:) ;
                    CIs_h2(4,:),    [NaN NaN],     CIs_e2(4,:) ;
                    CIs_h2(5,:),    [NaN NaN],     CIs_e2(5,:) ;
                    CIs_h2(6,:),    [NaN NaN],     CIs_e2(6,:) };
            % CIs = { CIs_h2(1,:),  [NaN NaN],   CIs_e2(1,:) ;
            %         CIs_h2(2,:),  [NaN NaN],   CIs_e2(2,:) ;
            %         Nc.CI_MZDZ*2, Nc.CI_DZSib, [NaN NaN]   ;
            %         Oc.CI_MZDZ*2, Oc.CI_DZSib, [NaN NaN]   ;
            %         CIs_h2(3,:),  [NaN NaN],   CIs_e2(3,:) ;
            %         CIs_h2(4,:),  [NaN NaN],   CIs_e2(4,:) ;
            %         CIs_h2(5,:),  [NaN NaN],   CIs_e2(5,:) ;
            %         CIs_h2(6,:),  [NaN NaN],   CIs_e2(6,:) };
            
    end
    
else
    
    CIs = repmat({NaN(1,2)},length(Names),3);

end


if ~ACEfit_Par.NoImg
    
    Names = {Names{:},'MaxLRT'};
    
    Ests = [ Ests;
             ACEfit_Par.mT, NaN, NaN ];
    
    CIs = [ CIs;
            {[NaN NaN],[NaN NaN],[NaN NaN]} ];
    
    if ACEfit_Par.nPerm>0 
        load(fullfile(ACEfit_Par.ResDir,'Pvals_Max_h2')) % FWE voxel-wise p-value p_T (and cluster-based p-value p_K & p_M)        
        Ps = [ Ps;
               p_T, NaN, NaN ];
    else
        Ps = [ Ps;
               NaN, NaN, NaN ];
    end
    
    if ACEfit_Par.Vs.ClustInf
        
        Names = {Names{:},'MaxSize','MaxMass'};
        
        Ests = [ Ests;
                 ACEfit_Par.mK, NaN, NaN;
                 ACEfit_Par.mM, NaN, NaN ];
        
        CIs = [ CIs;
                {[NaN NaN],[NaN NaN],[NaN NaN];
                 [NaN NaN],[NaN NaN],[NaN NaN]} ];
        
        if ACEfit_Par.nPerm>0            
            Ps = [ Ps;
                   p_K, NaN, NaN;
                   p_M, NaN, NaN ];  
        else
            Ps = [ Ps;
                   NaN, NaN, NaN;
                   NaN, NaN, NaN ];
        end
        
    end
    
end



if StdOut
    % Write to a string first so I can replace NaN's with empty space
    
    str = [...
           sprintf('%-10s    Heritability h^2 (aka a^2)               Common Env. c^2                  Error/Uniq. Env. e^2        \n','')...
           sprintf('%-10s----------------------------------  ----------------------------------  ----------------------------------\n','')...
           sprintf('          Estimate        95%%CI        P-val  Estimate        95%%CI        P-val  Estimate        95%%CI        P-val  \n')];
    for i=1:size(Ests,1)
        str = [str ...
               sprintf('%-10s %7.3f  ( %5.3f , %5.3f ) %6.3f   %7.3f  ( %5.3f , %5.3f ) %6.3f   %7.3f  ( %5.3f , %5.3f ) %6.3f\n',...
               Names{i},...
               Ests(i,1),CIs{i,1},Ps(i,1),...
               Ests(i,2),CIs{i,2},Ps(i,2),...
               Ests(i,3),CIs{i,3},Ps(i,3))];
    end
    
    str = strrep(str,' NaN','.   ');
    
    fprintf('%s',str);
    
else
    % Write to a string first so I can replace NaN's with empty cells
    
    str = [...
           sprintf(',Heritability h^2 (aka a^2),,,,,Common Env c^2,,,,,Error/Uniq. Env,,,\n')...
           sprintf(',Estimate,95%%CI lb,95%%CI ub,P-value,,Estimate,95%%CI lb,95%%CI ub,P-value,,Estimate,95%%CI lb,95%%CI ub,P-value\n')];
    for i=1:size(Ests,1)
        str = [str ...
               sprintf('%s,%f,%f,%f,%f,,%f,%f,%f,%f,,%f,%f,%f,%f\n',...
               Names{i},...
               Ests(i,1),CIs{i,1},Ps(i,1),...
               Ests(i,2),CIs{i,2},Ps(i,2),...
               Ests(i,3),CIs{i,3},Ps(i,3));];
    end
    
    str = strrep(str,'NaN','');
    
    fp = fopen(OutFile,'w');
    % fp = fopen([OutFile '.csv'],'w');
    fprintf(fp,'%s',str);
    fclose(fp);
    
end

if nargout>=1, varargout{1} = Ests; end
if nargout>=2, varargout{2} = CIs;  end
if nargout>=3, varargout{3} = Ps;   end
if nargout>=4, varargout{4} = {Names',{'a^2','c^2','e^2'}}; end

return

