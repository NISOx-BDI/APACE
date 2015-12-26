function [hh,HH] = LRSD_exact_ACE(nMZ,nDZ,n,Res,XTX,pinvXTX,sig2)
%
% Heritability estimation by LR-SD.
%

sdY    = zeros(3,1);
sdY(1) = sum((Res(    [1:2:nMZ])-Res(    [2:2:nMZ])).^2);
sdY(2) = sum((Res(nMZ+[1:2:nDZ])-Res(nMZ+[2:2:nDZ])).^2);
sdY(3) = (n^2-n)*sig2 - sdY(1) - sdY(2);

HH     = zeros(3,4);
    
XTy    = [ sdY(2)+2*sdY(3) ; 
           2*sdY(3)        ;
           2*sum(sdY)      ];
    
% XTy    = [ 0 1 2 
%            0 0 2
%            2 2 2 ]*y;

%
% HH = [E,AE,CE,ACE]
%
% hh = [0, 0, E]
HH(3,1)     = pinvXTX{1}*XTy(3);
% hh = [A, 0, E]
HH(1:2:3,2) = pinvXTX{2}*XTy(1:2:3);
% hh = [0, C, E]
HH(2:3,3)   = pinvXTX{3}*XTy(2:3);
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
