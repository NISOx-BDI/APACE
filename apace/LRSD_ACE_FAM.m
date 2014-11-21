function [hh,HH] = LRSD_ACE_FAM(I_MZ,I_DZ,I_Sib,Res,XTX,pinvXTX,sigma2)
%
% Heritability estimation by LR-SD using residuals
%
% INPUTS:
% I_MZ    - Index for MZ twin pairs
% I_DZ    - Index for DZ twin pairs
% I_Sib   - Index for sibling pairs
% Res     - Residuals from regressing the original data Y on X with OLS.
% XTX     - X'*X
% pinvXTX - pinv(X'*X)
% sigma2  - Residual sum of squares (RSS)
%
% OUTPUTS:
% hh      - Parameter estimate
% HH      - 4 possible estimates [E AE CE ACE]
%

% 4 possible models
MaxMod = 4;         
n      = length(Res);

% Data-dependent computations

ResMZ  = reshape(Res(I_MZ),size(I_MZ));
ResDZ  = reshape(Res(I_DZ),size(I_DZ));
ResSib = reshape(Res(I_Sib),size(I_Sib));

y    = zeros(3,1);
y(1) = sum((ResMZ(:,1)-ResMZ(:,2)).^2);
y(2) = sum((ResDZ(:,1)-ResDZ(:,2)).^2)+sum((ResSib(:,1)-ResSib(:,2)).^2);
y(3) = (n^2-n)*sigma2-y(1)-y(2);


HH  = zeros(3,MaxMod);

XTy = [ y(2)+2*y(3) 
        2*y(3) 
        2*(n^2-n)*sigma2 ];

% hh = [0, 0, E]
HH(3,1)     = pinvXTX{1}*XTy(3);

% hh = [A, 0, E]
est         = pinvXTX{2}*XTy(1:2:3);
HH(1:2:3,2) = est;

% hh = [0, C, E]
est         = pinvXTX{3}*XTy(2:3);
HH(2:3,3)   = est;

% hh = [A, C, E]
HH(:,4)     = pinvXTX{4}*XTy;


if min(HH(:,4))>=0
    hh  = HH(:,4);
else
    hhs = HH(:,min(HH)>=0);
    if ( min(HH(:,2))>=0 && min(HH(:,3))>=0 )
        fhh = diag(hhs'*XTX*hhs) - 2*hhs'*XTy;
        hhs = hhs(:,fhh==min(fhh));
    end
    hh  = hhs(:,end);
end

return
