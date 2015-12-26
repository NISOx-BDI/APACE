function PrepParallel(ACEfit_Par,nParallel)
%
% Update the input structure array
%
% nParallel - Number of parallel runs
%

if isempty(nParallel)
    nParallel = 1;
end

ACEfit_Par.RunID = 1:nParallel;

nPerm = ACEfit_Par.nPerm;
nMZF  = ACEfit_Par.nMZF;
nDZF  = ACEfit_Par.nDZF;

% Permutation preparation
if nPerm>0
    
    if rem(nPerm,nParallel)
        error('''nPerm'' cannot be divided evenly by ''nParallel''!!')
    end    
    ACEfit_Par.nPermPerRun = ones(nParallel,1)*nPerm/nParallel;
    
    ACEfit_Par.Perm_index  = CreatePerm(nMZF,nDZF,nPerm);
    
end

% Update and Save the main input file in the result folder
save(fullfile(ACEfit_Par.ResDir,'ACEfit_Par'),'ACEfit_Par','-v7.3');

return
