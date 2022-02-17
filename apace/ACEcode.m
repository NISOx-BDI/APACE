function x = ACEcode(ACE)
%
% Maps active (non-zero) ACE parameters to a integer code
%
%_______________________________________________________________________
% Version: http://github.com/NISOx-BDI/APACE/tree/$Format:%h$
%          $Format:%ci$

Code = [1 2 3];
x    = sum( Code(abs(ACE)>0) );

return

