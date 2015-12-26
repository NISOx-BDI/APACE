function [mT,uCount] = ACEfit_Perm(ACEfit_Par,Plabel)
%
% Permutations
%
% Perm_label - Permutation relabeling
%
nMZ    = ACEfit_Par.nMZF*2;
nDZ    = ACEfit_Par.nDZF*2;
X      = ACEfit_Par.X;
Y      = ACEfit_Par.Y;
I_data = ACEfit_Par.I_data;

load(fullfile(ACEfit_Par.ResDir,'Results'),'Sig2');

[nElm,n] = size(Y);

pinvXTX = X*pinv(X);
Pindex  = reshape([2*Plabel-1; 2*Plabel],nMZ+nDZ,1);

% Initalize result variable
[ERV,uCount] = deal(zeros(nElm));

for id = 1:length(I_data)
    
    i    = I_data(id); 
    
    iY   = Y(i,:)';
    iRes = iY - pinvXTX*iY;
    iStd = iRes/std(iRes);
    
    for j = I_data(id+1:end)
        
        jY   = Y(j,:)';
        jRes = jY - pinvXTX*jY;
        jStd = jRes/std(jRes);
        
        [ERV12] = Perm_ERV(nMZ,nDZ,n,iStd,jStd,1,1,reshape(ACEfit_Par.sdY(i,j,:),4,1),Pindex);
        % [ERV12] = Perm_ERV(nMZ,nDZ,n,iStd,jStd,var(iStd),var(jStd),reshape(ACEfit_Par.sdY(i,j,:),4,1),Pindex);
        
        ERV(i,j) = abs(ERV12);
        
    end
        
end

% Obtain the test statistic: ERV
Stats = ERV;
mT    = max(max(Stats));
uCount(Stats>=ACEfit_Par.Stats) = 1;

fprintf('.\n');

return

