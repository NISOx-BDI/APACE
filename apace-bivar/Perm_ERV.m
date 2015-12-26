function [ERV] = Perm_ERV(nMZ,nDZ,n,Res1,Res2,sig1,sig2,sdY,Pindex)
%
% Parameter estimation by LR-SD for each permutation.
%
nOTW = n^2-n-nMZ-nDZ;
N    = nMZ*nDZ+4*nMZ*nOTW+nDZ*nOTW;

Res1(1:nMZ+nDZ) = Res1(Pindex);
Res2(1:nMZ+nDZ) = Res2(Pindex);

sdY(2) = sum((Res1(    [1:2:nMZ])-Res2(    [2:2:nMZ])).^2) + sum((Res1(    [2:2:nMZ])-Res2(    [1:2:nMZ])).^2);
sdY(3) = sum((Res1(nMZ+[1:2:nDZ])-Res2(nMZ+[2:2:nDZ])).^2) + sum((Res1(nMZ+[2:2:nDZ])-Res2(nMZ+[1:2:nDZ])).^2);

ERV  = (-(nDZ+2*nOTW)*sdY(2)+(nMZ-nOTW)*sdY(3)+(2*nMZ+nDZ)*sdY(4))/N/sqrt(sig1*sig2); 

ERV  = min(max(ERV,-1),1);

return
