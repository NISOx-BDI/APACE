function x = ACEcode(ACE)
%
% Maps active (non-zero) ACE parameters to a integer code
%

Code = [1 2 3];
x    = sum( Code(abs(ACE)>0) );

return

