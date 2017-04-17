function [pinvZTZ] = Mk_ZTZpinv(nMZSD,nDZSD,nSibSD,nOTW,n)
%
% Compute ZTZ = Z'*Z & inverse matrix of ZTZ
%
% Z = [ 0 0 2
%        ...
%       1 0 2
%        ...
%       2 2 2
%        ...  ]
%
% ZTZ = [     nDZSD+nSibSD+4*nOTW 4*nOTW 2*(nDZSD+nSibSD)+4*nOTW
%                          4*nOTW 4*nOTW                  4*nOTW
%         2*(nDZSD+nSibSD)+4*nOTW 4*nOTW               2*(n^2-n) ];
%
%_______________________________________________________________________
% Version: http://github.com/nicholst/APACE/tree/$Format:%h$
%          $Format:%ci$

pinvZTZ    = cell(4,1);

%pinvZTZ{1} = 1/2/(n^2-n);

fixterm    = 2*nMZSD*(2*nDZSD+2*nSibSD+8*nOTW)+4*(nDZSD+nSibSD)*nOTW;
pinvZTZ{2} = [                2*(n^2-n) -2*nDZSD-2*nSibSD-4*nOTW            
               -2*nDZSD-2*nSibSD-4*nOTW      nDZSD+nSibSD+4*nOTW ]/fixterm;

fixterm    = 1/4/(nMZSD+nDZSD+nSibSD);
pinvZTZ{3} = [ 1/4/nOTW+fixterm -fixterm
                       -fixterm  fixterm ];

diagACE    = diag([1/nMZSD+1/(nDZSD+nSibSD) 1/4/nOTW+1/4/nMZSD+1/(nDZSD+nSibSD) 1/4/nMZSD]);
triuACE    = [ 0 -1/2/nMZSD-1/(nDZSD+nSibSD) -1/2/nMZSD
               0                           0  1/4/nMZSD
               0                           0          0 ];
pinvZTZ{4} = diagACE+triuACE+triuACE';


% pinvZTZ = { 1/ZTZ(3,3),...
%             pinv(ZTZ(1:2:3,1:2:3)),...
%             pinv(ZTZ(2:3,2:3)),...
%             pinv(ZTZ) };

return
