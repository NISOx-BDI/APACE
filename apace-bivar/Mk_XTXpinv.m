function [XTX,pinvXTX] = Mk_XTXpinv(nMZSD,nDZSD,nSibSD,n)
%
% Compute XTX = X'*X & inverse matrix of XTX
%
% X = [ 0 0 2
%        ...
%       1 0 2
%        ...
%       2 2 2
%        ...  ]
%

nOTW = (n^2-n)/2-nMZSD-nDZSD-nSibSD;

XTX = [ nDZSD+nSibSD+4*nOTW     4*nOTW 2*(nDZSD+nSibSD)+4*nOTW
        4*nOTW                  4*nOTW 4*nOTW
        2*(nDZSD+nSibSD)+4*nOTW 4*nOTW 2*(n^2-n)               ];

pinvXTX    = cell(4,1);

pinvXTX{1} = 1/XTX(3,3);

deno       = XTX(1,3)^2 - XTX(1,1)*XTX(3,3);
pinvXTX{2} = [ -XTX(3,3)  XTX(1,3)
                XTX(1,3) -XTX(1,1) ]/deno;

deno       = XTX(2,3)^2 - XTX(2,2)*XTX(3,3);
pinvXTX{3} = [ -XTX(3,3)  XTX(2,3)
                XTX(2,3) -XTX(2,2) ]/deno;

deno       = XTX(1,2)*XTX(3,3)-2*XTX(1,2)*XTX(1,3)+XTX(1,3)^2-XTX(1,1)*XTX(3,3)+XTX(1,1)*XTX(1,2);
diagACE    = [ XTX(1,2)-XTX(3,3) 0                                       0 
               0                 (XTX(1,3)^2-XTX(1,1)*XTX(3,3))/XTX(1,2) 0
               0                 0                                       XTX(1,2)-XTX(1,1) ];
triuACE    = [ 0 XTX(3,3)-XTX(1,3) XTX(1,3)-XTX(1,2)
               0 0                 XTX(1,1)-XTX(1,3)
               0 0                 0                 ];
pinvXTX{4} = (diagACE+triuACE+triuACE')/deno;


% pinvXTX = { 1/XTX(3,3),...
%             pinv(XTX(1:2:3,1:2:3)),...
%             pinv(XTX(2:3,2:3)),...
%             pinv(XTX) };

return
