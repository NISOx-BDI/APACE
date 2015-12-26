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

nElm = size(Y,1);

% Make residuals, compute total variance, initialise
Res           = Y - Y*(X*pinv(X))';
Sigma2        = sum(Res.^2,2)/(nB-size(X,2));
[XTX,pinvXTX] = Mk_XTXpinv(size(I_MZ,1),size(I_DZ,1),size(I_Sib,1),nB);


switch upper(ACEfit_Par.Model)
    
    case 'ACE'
        
        % Initalize results variables
        [h2,c2,e2] = deal(zeros(nElm,1));
        
        for i = I_data
            
            iRes = Res(i,:)';
            sig2 = Sigma2(i);
            
            %
            % Heritability Estimation
            %
            [hha,~] = LRSD_ACE_FAM(I_MZ,I_DZ,I_Sib,iRes,XTX,pinvXTX,sig2);
            hha     = hha/sum(hha);
            
            h2(i)   = hha(1);
            c2(i)   = hha(2);
            e2(i)   = hha(3);
            
        end
        
        h2       = h2(I_data);
        c2       = c2(I_data);
        e2       = e2(I_data);
        Sigma2   = Sigma2(I_data);
        
        Aboot    = Sigma2.*h2;
        wh2      = mean(Aboot)/mean(Sigma2);
        
        quartile = quantile(h2,[0.50 0.75]);
        medh2    = quartile(1);
        q3h2     = quartile(2);
        meanh2   = mean(h2);
        mGmedh2  = mean(h2(h2>medh2));
        % if isnan(mGmedh2)
        %     mGmedh2  = mean(h2(h2>=medh2));
        % end
        mGq3h2   = mean(h2(h2>q3h2));
        % if isnan(mGq3h2)
        %     mGq3h2   = mean(h2(h2>=q3h2));
        % end
        
        % Summary statistics: meanh2, wh2, median, q3, mean(h2>median), mean(h2>q3)
        SummaryA = [meanh2; wh2; medh2; q3h2; mGmedh2; mGq3h2];
        
        Cboot    = Sigma2.*c2;
        wc2      = mean(Cboot)/mean(Sigma2);
        
        quartile = quantile(c2,[0.50 0.75]);
        medc2    = quartile(1);
        q3c2     = quartile(2);
        meanc2   = mean(c2);
        mGmedc2  = mean(c2(c2>medc2));
        % if isnan(mGmedc2)
        %     mGmedc2  = mean(c2(c2>=medc2));
        % end
        mGq3c2   = mean(c2(c2>q3c2));
        % if isnan(mGq3c2)
        %     mGq3c2   = mean(c2(c2>=q3c2));
        % end
        
        % Summary statistics: meanc2, wc2, median, q3, mean(c2>median), mean(c2>q3)
        SummaryC = [meanc2; wc2; medc2; q3c2; mGmedc2; mGq3c2];
        
        Eboot    = Sigma2.*e2;
        we2      = mean(Eboot)/mean(Sigma2);
        
        quartile = quantile(e2,[0.50 0.75]);
        mede2    = quartile(1);
        q3e2     = quartile(2);
        meane2   = mean(e2);
        mGmede2  = mean(e2(e2>mede2));
        % if isnan(mGmede2)
        %     mGmede2  = mean(e2(e2>=mede2));
        % end
        mGq3e2   = mean(e2(e2>q3e2));
        % if isnan(mGq3e2)
        %     mGq3e2   = mean(e2(e2>=q3e2));
        % end
        
        % Summary statistics: meane2, we2, median, q3, mean(e2>median), mean(e2>q3)
        SummaryE = [meane2; we2; mede2; q3e2; mGmede2; mGq3e2];
        
        % Summary statistics for h2, c2 and e2
        SummaryACE = [SummaryA SummaryC SummaryE];
        
    case 'AE'
                
        % Initalize results variables
        [h2,e2] = deal(zeros(nElm,1));
        
        for i = I_data
            
            iRes = Res(i,:)';
            sig2 = Sigma2(i);
            
            %
            % Heritability Estimation
            %
            [~,HH] = LRSD_ACE_FAM(I_MZ,I_DZ,I_Sib,iRes,XTX,pinvXTX,sig2);
            if ( min(HH(:,2))>=0 )
                hha = HH(:,2);
            else
                hha = HH(:,1);
            end
            hha     = hha/sum(hha);
            
            h2(i)   = hha(1);
            e2(i)   = hha(3);
            
        end
        
        h2       = h2(I_data);
        e2       = e2(I_data);
        Sigma2   = Sigma2(I_data);
        
        Aboot    = Sigma2.*h2;
        wh2      = mean(Aboot)/mean(Sigma2);
        
        quartile = quantile(h2,[0.50 0.75]);
        medh2    = quartile(1);
        q3h2     = quartile(2);
        meanh2   = mean(h2);
        mGmedh2  = mean(h2(h2>medh2));
        % if isnan(mGmedh2)
        %     mGmedh2  = mean(h2(h2>=medh2));
        % end
        mGq3h2   = mean(h2(h2>q3h2));
        % if isnan(mGq3h2)
        %     mGq3h2   = mean(h2(h2>=q3h2));
        % end
        
        % Summary statistics: meanh2, wh2, median, q3, mean(h2>median), mean(h2>q3)
        SummaryA = [meanh2; wh2; medh2; q3h2; mGmedh2; mGq3h2];
        
        Eboot    = Sigma2.*e2;
        we2      = mean(Eboot)/mean(Sigma2);
        
        quartile = quantile(e2,[0.50 0.75]);
        mede2    = quartile(1);
        q3e2     = quartile(2);
        meane2   = mean(e2);
        mGmede2  = mean(e2(e2>mede2));
        % if isnan(mGmede2)
        %     mGmede2  = mean(e2(e2>=mede2));
        % end
        mGq3e2   = mean(e2(e2>q3e2));
        % if isnan(mGq3e2)
        %     mGq3e2   = mean(e2(e2>=q3e2));
        % end
        
        % Summary statistics: meane2, we2, median, q3, mean(e2>median), mean(e2>q3)
        SummaryE = [meane2; we2; mede2; q3e2; mGmede2; mGq3e2];
        
        % Summary statistics for h2, c2 and e2
        SummaryACE = [SummaryA SummaryE];
        
end

return
