function [hh,hhAE,hhCE] = LRSD_ACE_VEC(I_MZ,I_DZ,I_Sib,Res,sigma2)
%
% Heritability estimation using LR-SD with residuals
%
% INPUTS:
% I_MZ    - Index for MZ twin pairs
% I_DZ    - Index for DZ twin pairs
% I_Sib   - Index for sibling pairs
% Res     - Residuals from regressing the original data Y on X with OLS
% sigma2  - variance estimator associated with residual sum of squares (RSS)
%
% OUTPUTS:
% hh      - Variance parameter estimate
% hh*     - The parameter estimates of the possible models [E AE CE ACE]
%
%_______________________________________________________________________
% Version: http://github.com/NISOx-BDI/APACE/tree/$Format:%h$
%          $Format:%ci$

[nElm,n]  = size(Res);
nMZSD     = size(I_MZ,1);
nDZSD     = size(I_DZ,1);
nSibSD    = size(I_Sib,1);
nOTW      = (n^2-n)/2-nMZSD-nDZSD-nSibSD;
[pinvZTZ] = Mk_ZTZpinv(nMZSD,nDZSD,nSibSD,nOTW,n);

% Data-dependent computations

ResMZ1  = Res(:,I_MZ(:,1));
ResDZ1  = Res(:,I_DZ(:,1));
ResSib1 = Res(:,I_Sib(:,1));

ResMZ2  = Res(:,I_MZ(:,2));
ResDZ2  = Res(:,I_DZ(:,2));
ResSib2 = Res(:,I_Sib(:,2));

y      = zeros(nElm,3);
y(:,1) = sum((ResMZ1-ResMZ2).^2,2);
y(:,2) = sum((ResDZ1-ResDZ2).^2,2)+sum((ResSib1-ResSib2).^2,2);
y(:,3) = (n^2-n)*sigma2-y(:,1)-y(:,2);


[hh,hhAE,hhCE] = deal(repmat([0 0 1],nElm,1));

ZTD = [ y(:,2)+2*y(:,3) 2*y(:,3) 2*(n^2-n)*sigma2 ];

% hhE   -- [0, 0, E]
% hhE(:,3)      = XTy(:,3)*pinvZTZ{1};

% hhAE  -- [A, 0, E]
hhAE(:,1:2:3) = ZTD(:,1:2:3)*pinvZTZ{2};

% hhCE  -- [0, C, E]
hhCE(:,2:3)   = ZTD(:,2:3)*pinvZTZ{3};

% hhACE -- [A, C, E]
hhACE         = ZTD*pinvZTZ{4};


hhACEplus  = min(hhACE,[],2)>=0;
hhACEminus = min(hhACE,[],2)<0;
hhAEplus   = min(hhAE,[],2)>=0;
hhAEminus  = min(hhAE,[],2)<0;
hhCEplus   = min(hhCE,[],2)>=0;
hhCEminus  = min(hhCE,[],2)<0;

hh(hhACEplus,:) = hhACE(hhACEplus,:);
hhplus          = hhACEminus & hhAEplus & hhCEminus;
hh(hhplus,:)    = hhAE(hhplus,:);
hhplus          = hhACEminus & hhAEminus & hhCEplus;
hh(hhplus,:)    = hhCE(hhplus,:);

hhplus             = hhACEminus & hhAEplus & hhCEplus;
hhAEs              = hhAE(hhplus,:);
hhCEs              = hhCE(hhplus,:);
ZTDs               = ZTD(hhplus,:);
RssAE              = (nDZSD+nSibSD+4*nOTW)*(hhAEs(:,1).^2) + ...
                     (4*nDZSD+4*nSibSD+8*nOTW)*(hhAEs(:,1).*hhAEs(:,3)) + ...
                     2*(n^2-n)*(hhAEs(:,3).^2) - ...
                     2*(hhAEs(:,1).*ZTDs(:,1)) - ...
                     2*(hhAEs(:,3).*ZTDs(:,3));
RssCE              = 4*nOTW*(hhCEs(:,2).^2) + ...
                     8*nOTW*(hhCEs(:,2).*hhCEs(:,3)) + ...
                     2*(n^2-n)*(hhCEs(:,3).^2) - ...
                     2*(hhCEs(:,2).*ZTDs(:,2)) - ...
                     2*(hhCEs(:,3).*ZTDs(:,3));
hhs                = hhAEs;
hhs(RssAE>RssCE,:) = hhCEs(RssAE>RssCE,:);
hh(hhplus,:)       = hhs;

hhAE(hhAEminus,1:2) = 0;
hhAE(hhAEminus,3)   = 1;

hhCE(hhCEminus,1:2) = 0;
hhCE(hhCEminus,3)   = 1;


return
