function [sdY,rhoP,ERV] = LRSD_BIVAR(nMZ,nDZ,n,Res1,Res2,sig1,sig2)
%
% Parameter estimation by LR-SD.
%

nOTW = n^2-n-nMZ-nDZ;
N1   = n*nOTW;
N2   = nMZ*nDZ+4*nMZ*nOTW+nDZ*nOTW;
% N2   = nMZ*nDZ+nMZ*nOTW+nDZ*nOTW;

sdY    = zeros(4,1);
sdY(1) = sum((Res1-Res2).^2);
sdY(2) = sum((Res1(    [1:2:nMZ])-Res2(    [2:2:nMZ])).^2) + sum((Res1(    [2:2:nMZ])-Res2(    [1:2:nMZ])).^2);
sdY(3) = sum((Res1(nMZ+[1:2:nDZ])-Res2(nMZ+[2:2:nDZ])).^2) + sum((Res1(nMZ+[2:2:nDZ])-Res2(nMZ+[1:2:nDZ])).^2);

sum12 = 0;
for i = 1:n
    sum12 = sum12 + sum((Res1(i)-Res2).^2);
end
sdY(4) = sum12 - sdY(1) - sdY(2) - sdY(3);

Par1 = (n*sdY(4)-nOTW*sdY(1))/(2*N1);
Par2 = (-(nDZ+2*nOTW)*sdY(2)+(nMZ-nOTW)*sdY(3)+(2*nMZ+nDZ)*sdY(4))/N2;
% Par1 = (n*y(4)-nOTW*y(1))/(2*N1);
% Par2 = (-(nDZ+nOTW/2)*y(2)+(nMZ-nOTW)*y(3)+(nMZ/2+nDZ)*y(4))/N2;

rhoP = Par1/sqrt(sig1*sig2);
ERV  = Par2/sqrt(sig1*sig2); 

rhoP = min(max(rhoP,-1),1);
ERV  = min(max(ERV,-1),1);

return
