function [B_CI_LB,B_CI_UB] = Bootstrap_CI(ST,alpha)
%
% Generate Bootstrap CI.
%
% ST     - bootstrap estimates and the original estimate
% alpha  - level
%
%_______________________________________________________________________
% Version: http://github.com/NISOx-BDI/APACE/tree/$Format:%h$
%          $Format:%ci$

% 1-alpha/2 quantile of the standard normal distribution
ST_alpha = spm_invNcdf(1-alpha/2);

nBoot    = length(ST)-1;
STO      = ST(end);       % Original estimate
ST       = ST(1:nBoot);   % Bootstrap estimates

% 1) Standard CI
stdB     = std(ST);
CI1L     = STO - ST_alpha*stdB;
CI1U     = STO + ST_alpha*stdB;

% 2) Percentile CI
sSTB     = sort(ST);
CI2L     = sSTB(min(max(1, round((nBoot+1)*alpha/2)     ),end));
CI2U     = sSTB(min(max(1, round((nBoot+1)*(1-alpha/2)) ),end));

% 3) Bias corrected percentile CI
BCa_z    = spm_invNcdf(sum(ST<=STO)/(nBoot+1));
BC_aL    = spm_Ncdf(2*BCa_z - ST_alpha);
BC_aU    = spm_Ncdf(2*BCa_z + ST_alpha);
CI3L     = sSTB(min(max(1,round(nBoot*BC_aL)),end));
CI3U     = sSTB(min(max(1,round(nBoot*BC_aU)),end));

B_CI_LB  = min([CI1L, CI2L, CI3L]);
B_CI_UB  = max([CI1U, CI2U, CI3U]);

return
