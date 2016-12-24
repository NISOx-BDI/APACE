function [SummaryACE] = ACEfit_Boot(ACEfit_Par,Boot_label,Boot_kin)
%
% Bootstrapping.
%

nB       = size(Boot_kin,1);
Boot_kin = triu(Boot_kin);
fk_MZ    = find(Boot_kin==1);
fk_DZ    = find(Boot_kin==2);
fk_Sib   = find(Boot_kin==3);
I_MZ     = [rem(fk_MZ,nB)  (fk_MZ  - rem(fk_MZ,nB) )/nB+1];
I_DZ     = [rem(fk_DZ,nB)  (fk_DZ  - rem(fk_DZ,nB) )/nB+1];
I_Sib    = [rem(fk_Sib,nB) (fk_Sib - rem(fk_Sib,nB))/nB+1];
I_data   = ACEfit_Par.I_data;

Y = ACEfit_Par.Y;
Y = Y(:,Boot_label);

X = ACEfit_Par.X;
X = X(Boot_label,:);

% Make residuals, compute total variance, initialise
Res    = Y - Y*(X*pinv(X))';
Res    = Res(I_data,:);
Sigma2 = sum(Res.^2,2)/(nB-size(X,2));

%
% Heritability estimation
%
[hhACE0,hhAE0,~] = LRSD_ACE_VEC(I_MZ,I_DZ,I_Sib,Res,Sigma2);

% Hack, for possibility of all-constant Bootstrap sample with discrete data:
% Set ACE to [0 0 1]
Iconst = find(Sigma2==0);
hhACE0(Iconst,1:2)=0; hhACE0(Iconst,3)=1;
 hhAE0(Iconst,1:2)=0;  hhAE0(Iconst,3)=1;

switch upper(ACEfit_Par.Model)
    
    case 'ACE'
        
        hh       = hhACE0./repmat(sum(hhACE0,2),1,3);
        
        h2       = hh(:,1);
        c2       = hh(:,2);
        e2       = hh(:,3);
        
        Aboot    = Sigma2.*h2;
        wh2      = mean(Aboot)/mean(Sigma2);
        
        quartile = quantile(h2,[0.50 0.75]);
        medh2    = quartile(1);
        q3h2     = quartile(2);
        meanh2   = mean(h2);
        mGmedh2  = mean(h2(h2>=medh2));
        mGq3h2   = mean(h2(h2>=q3h2));
        
        % Summary statistics: meanh2, wh2, median, q3, mean(h2>median), mean(h2>q3)
        SummaryA = [meanh2; wh2; medh2; q3h2; mGmedh2; mGq3h2];
        
        Cboot    = Sigma2.*c2;
        wc2      = mean(Cboot)/mean(Sigma2);
        
        quartile = quantile(c2,[0.50 0.75]);
        medc2    = quartile(1);
        q3c2     = quartile(2);
        meanc2   = mean(c2);
        mGmedc2  = mean(c2(c2>=medc2));
        mGq3c2   = mean(c2(c2>=q3c2));
        
        % Summary statistics: meanc2, wc2, median, q3, mean(c2>median), mean(c2>q3)
        SummaryC = [meanc2; wc2; medc2; q3c2; mGmedc2; mGq3c2];
        
        Eboot    = Sigma2.*e2;
        we2      = mean(Eboot)/mean(Sigma2);
        
        quartile = quantile(e2,[0.50 0.75]);
        mede2    = quartile(1);
        q3e2     = quartile(2);
        meane2   = mean(e2);
        mGmede2  = mean(e2(e2>=mede2));
        mGq3e2   = mean(e2(e2>=q3e2));
        
        % Summary statistics: meane2, we2, median, q3, mean(e2>median), mean(e2>q3)
        SummaryE = [meane2; we2; mede2; q3e2; mGmede2; mGq3e2];
        
        % Summary statistics for h2, c2 and e2
        SummaryACE = [SummaryA SummaryC SummaryE];
        
    case 'AE'
        
        hh       = hhAE0./repmat(sum(hhAE0,2),1,3);
        
        h2       = hh(:,1);
        e2       = hh(:,3);
        
        Aboot    = Sigma2.*h2;
        wh2      = mean(Aboot)/mean(Sigma2);
        
        quartile = quantile(h2,[0.50 0.75]);
        medh2    = quartile(1);
        q3h2     = quartile(2);
        meanh2   = mean(h2);
        mGmedh2  = mean(h2(h2>=medh2));
        mGq3h2   = mean(h2(h2>=q3h2));
        
        % Summary statistics: meanh2, wh2, median, q3, mean(h2>median), mean(h2>q3)
        SummaryA = [meanh2; wh2; medh2; q3h2; mGmedh2; mGq3h2];
        
        Eboot    = Sigma2.*e2;
        we2      = mean(Eboot)/mean(Sigma2);
        
        quartile = quantile(e2,[0.50 0.75]);
        mede2    = quartile(1);
        q3e2     = quartile(2);
        meane2   = mean(e2);
        mGmede2  = mean(e2(e2>=mede2));
        mGq3e2   = mean(e2(e2>=q3e2));
        
        % Summary statistics: meane2, we2, median, q3, mean(e2>median), mean(e2>q3)
        SummaryE = [meane2; we2; mede2; q3e2; mGmede2; mGq3e2];
        
        % Summary statistics for h2 and e2
        SummaryACE = [SummaryA SummaryE];
        
end

return
