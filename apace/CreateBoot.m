function Boot_index = CreateBoot(nG1,nG2,nG3,nBoot)
%
% Create balanced bootstrapping resampling labels
%
% nG1    - number of subjects or subject pairs or families in group 1
% nG2    - number of subjects or subject pairs or families in group 2
% nG3    - number of subjects or subject pairs or families in group 3
% nBoot  - number of bootstrapping replicates
%

% Derive seed from time
rng('shuffle');

% Balanced bootstrap
label_G1 = repmat(1:nG1, 1, nBoot)';
label_G2 = repmat(1:nG2, 1, nBoot)';
label_G3 = repmat(1:nG3, 1, nBoot)';
G1boot   = reshape(label_G1(randperm(nG1*nBoot)), nG1, nBoot)';
G2boot   = reshape(label_G2(randperm(nG2*nBoot)), nG2, nBoot)';
G3boot   = reshape(label_G3(randperm(nG3*nBoot)), nG3, nBoot)';

Boot_index = [G1boot nG1+G2boot nG1+nG2+G3boot];

return