function Boot_CIs(ACEfit_Par,varargin)
% FORMAT Boot_CIs(ACEfit_Par[,Balpha])
%
% Generate Bootstrap CIs.
%

if nargin>=2
    Balpha = varargin{1};
else
    Balpha = 0.05;
end

load(fullfile(ACEfit_Par.ResDir,'ACEfit_Boot'));

switch upper(ACEfit_Par.Model)
    
    case 'ACE'
        
        [CIs_h2,CIs_c2,CIs_e2] = deal(zeros(6,2));
        
        %
        % (1) mean, wh2, median, q3, mean(h2>median), mean(h2>q3)
        %
        [CIs_h2(1,1), CIs_h2(1,2)] = Bootstrap_CI(meanh2_ACE,  Balpha);
        [CIs_h2(2,1), CIs_h2(2,2)] = Bootstrap_CI(wh2_ACE,     Balpha);
        [CIs_h2(3,1), CIs_h2(3,2)] = Bootstrap_CI(medh2_ACE,   Balpha);
        [CIs_h2(4,1), CIs_h2(4,2)] = Bootstrap_CI(q3h2_ACE,    Balpha);
        [CIs_h2(5,1), CIs_h2(5,2)] = Bootstrap_CI(mGmedh2_ACE, Balpha);
        [CIs_h2(6,1), CIs_h2(6,2)] = Bootstrap_CI(mGq3h2_ACE,  Balpha);
        
        %
        % (2) mean, wc2, median, q3, mean(c2>median), mean(c2>q3)
        %
        [CIs_c2(1,1), CIs_c2(1,2)] = Bootstrap_CI(meanc2_ACE,  Balpha);
        [CIs_c2(2,1), CIs_c2(2,2)] = Bootstrap_CI(wc2_ACE,     Balpha);
        [CIs_c2(3,1), CIs_c2(3,2)] = Bootstrap_CI(medc2_ACE,   Balpha);
        [CIs_c2(4,1), CIs_c2(4,2)] = Bootstrap_CI(q3c2_ACE,    Balpha);
        [CIs_c2(5,1), CIs_c2(5,2)] = Bootstrap_CI(mGmedc2_ACE, Balpha);
        [CIs_c2(6,1), CIs_c2(6,2)] = Bootstrap_CI(mGq3c2_ACE,  Balpha);
        
        %
        % (3) mean, we2, median, q3, mean(e2>median), mean(e2>q3)
        %
        [CIs_e2(1,1), CIs_e2(1,2)] = Bootstrap_CI(meane2_ACE,  Balpha);
        [CIs_e2(2,1), CIs_e2(2,2)] = Bootstrap_CI(we2_ACE,     Balpha);
        [CIs_e2(3,1), CIs_e2(3,2)] = Bootstrap_CI(mede2_ACE,   Balpha);
        [CIs_e2(4,1), CIs_e2(4,2)] = Bootstrap_CI(q3e2_ACE,    Balpha);
        [CIs_e2(5,1), CIs_e2(5,2)] = Bootstrap_CI(mGmede2_ACE, Balpha);
        [CIs_e2(6,1), CIs_e2(6,2)] = Bootstrap_CI(mGq3e2_ACE,  Balpha);
        
        CIs_h2 = max(0,CIs_h2);
        CIs_h2 = min(1,CIs_h2);
        
        CIs_c2 = max(0,CIs_c2);
        CIs_c2 = min(1,CIs_c2);
        
        CIs_e2 = max(0,CIs_e2);
        CIs_e2 = min(1,CIs_e2);
        
        fprintf('The %.1f%% confidence interval for mean of h2 is (%.2f, %.2f). \n',                    100*(1-Balpha), CIs_h2(1,:));
        fprintf('The %.1f%% confidence interval for weighted mean of h2 is (%.2f, %.2f). \n',           100*(1-Balpha), CIs_h2(2,:));
        fprintf('The %.1f%% confidence interval for median (Q2) of h2 is (%.2f, %.2f). \n',             100*(1-Balpha), CIs_h2(3,:));
        fprintf('The %.1f%% confidence interval for the third quantile (Q3) of h2 is (%.2f, %.2f). \n', 100*(1-Balpha), CIs_h2(4,:));
        fprintf('The %.1f%% confidence interval for mean of h2>Q2 is (%.2f, %.2f). \n',                 100*(1-Balpha), CIs_h2(5,:));
        fprintf('The %.1f%% confidence interval for mean of h2>Q3 is (%.2f, %.2f). \n \n',              100*(1-Balpha), CIs_h2(6,:));
        
        fprintf('The %.1f%% confidence interval for mean of c2 is (%.2f, %.2f). \n',                    100*(1-Balpha), CIs_c2(1,:));
        fprintf('The %.1f%% confidence interval for weighted mean of c2 is (%.2f, %.2f). \n',           100*(1-Balpha), CIs_c2(2,:));
        fprintf('The %.1f%% confidence interval for median (Q2) of c2 is (%.2f, %.2f). \n',             100*(1-Balpha), CIs_c2(3,:));
        fprintf('The %.1f%% confidence interval for the third quantile (Q3) of c2 is (%.2f, %.2f). \n', 100*(1-Balpha), CIs_c2(4,:));
        fprintf('The %.1f%% confidence interval for mean of c2>Q2 is (%.2f, %.2f). \n',                 100*(1-Balpha), CIs_c2(5,:));
        fprintf('The %.1f%% confidence interval for mean of c2>Q3 is (%.2f, %.2f). \n \n',              100*(1-Balpha), CIs_c2(6,:));
        
        fprintf('The %.1f%% confidence interval for mean of e2 is (%.2f, %.2f). \n',                    100*(1-Balpha), CIs_e2(1,:));
        fprintf('The %.1f%% confidence interval for weighted mean of e2 is (%.2f, %.2f). \n',           100*(1-Balpha), CIs_e2(2,:));
        fprintf('The %.1f%% confidence interval for median (Q2) of e2 is (%.2f, %.2f). \n',             100*(1-Balpha), CIs_e2(3,:));
        fprintf('The %.1f%% confidence interval for the third quantile (Q3) of e2 is (%.2f, %.2f). \n', 100*(1-Balpha), CIs_e2(4,:));
        fprintf('The %.1f%% confidence interval for mean of e2>Q2 is (%.2f, %.2f). \n',                 100*(1-Balpha), CIs_e2(5,:));
        fprintf('The %.1f%% confidence interval for mean of e2>Q3 is (%.2f, %.2f). \n',                 100*(1-Balpha), CIs_e2(6,:));
        
        save(fullfile(ACEfit_Par.ResDir,'Boot_CIs'),'Balpha','CIs_h2','CIs_c2','CIs_e2');
        
    case 'AE'
        
        [CIs_h2,CIs_e2] = deal(zeros(6,2));
        
        %
        % (1) meanh2, wh2, median, q3, mean(h2>median), mean(h2>q3)
        %
        [CIs_h2(1,1), CIs_h2(1,2)] = Bootstrap_CI(meanh2_ACE,  Balpha);
        [CIs_h2(2,1), CIs_h2(2,2)] = Bootstrap_CI(wh2_ACE,     Balpha);
        [CIs_h2(3,1), CIs_h2(3,2)] = Bootstrap_CI(medh2_ACE,   Balpha);
        [CIs_h2(4,1), CIs_h2(4,2)] = Bootstrap_CI(q3h2_ACE,    Balpha);
        [CIs_h2(5,1), CIs_h2(5,2)] = Bootstrap_CI(mGmedh2_ACE, Balpha);
        [CIs_h2(6,1), CIs_h2(6,2)] = Bootstrap_CI(mGq3h2_ACE,  Balpha);
        
        %
        % (2) meane2, we2, median, q3, mean(e2>median), mean(e2>q3)
        %
        [CIs_e2(1,1), CIs_e2(1,2)] = Bootstrap_CI(meane2_ACE,  Balpha);
        [CIs_e2(2,1), CIs_e2(2,2)] = Bootstrap_CI(we2_ACE,     Balpha);
        [CIs_e2(3,1), CIs_e2(3,2)] = Bootstrap_CI(mede2_ACE,   Balpha);
        [CIs_e2(4,1), CIs_e2(4,2)] = Bootstrap_CI(q3e2_ACE,    Balpha);
        [CIs_e2(5,1), CIs_e2(5,2)] = Bootstrap_CI(mGmede2_ACE, Balpha);
        [CIs_e2(6,1), CIs_e2(6,2)] = Bootstrap_CI(mGq3e2_ACE,  Balpha);
        
        CIs_h2 = max(0,CIs_h2);
        CIs_h2 = min(1,CIs_h2);
        
        CIs_e2 = max(0,CIs_e2);
        CIs_e2 = min(1,CIs_e2);
        
        fprintf('The %.1f%% confidence interval for mean of h2 is (%.2f, %.2f). \n',                    100*(1-Balpha), CIs_h2(1,:));
        fprintf('The %.1f%% confidence interval for weighted mean of h2 is (%.2f, %.2f). \n',           100*(1-Balpha), CIs_h2(2,:));
        fprintf('The %.1f%% confidence interval for median (Q2) of h2 is (%.2f, %.2f). \n',             100*(1-Balpha), CIs_h2(3,:));
        fprintf('The %.1f%% confidence interval for the third quantile (Q3) of h2 is (%.2f, %.2f). \n', 100*(1-Balpha), CIs_h2(4,:));
        fprintf('The %.1f%% confidence interval for mean of h2>Q2 is (%.2f, %.2f). \n',                 100*(1-Balpha), CIs_h2(5,:));
        fprintf('The %.1f%% confidence interval for mean of h2>Q3 is (%.2f, %.2f). \n \n',              100*(1-Balpha), CIs_h2(6,:));
        
        fprintf('The %.1f%% confidence interval for mean of e2 is (%.2f, %.2f). \n',                    100*(1-Balpha), CIs_e2(1,:));
        fprintf('The %.1f%% confidence interval for weighted mean of e2 is (%.2f, %.2f). \n',           100*(1-Balpha), CIs_e2(2,:));
        fprintf('The %.1f%% confidence interval for median (Q2) of e2 is (%.2f, %.2f). \n',             100*(1-Balpha), CIs_e2(3,:));
        fprintf('The %.1f%% confidence interval for the third quantile (Q3) of e2 is (%.2f, %.2f). \n', 100*(1-Balpha), CIs_e2(4,:));
        fprintf('The %.1f%% confidence interval for mean of e2>Q2 is (%.2f, %.2f). \n',                 100*(1-Balpha), CIs_e2(5,:));
        fprintf('The %.1f%% confidence interval for mean of e2>Q3 is (%.2f, %.2f). \n',                 100*(1-Balpha), CIs_e2(6,:));
        
        save(fullfile(ACEfit_Par.ResDir,'Boot_CIs'),'Balpha','CIs_h2','CIs_e2');
        
end

return



