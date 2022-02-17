function PrepParallel(ACEfit_Par)
%
% Update the input structure array 'ACEfit_Par' for parallelization
%
%_______________________________________________________________________
% Version: http://github.com/NISOx-BDI/APACE/tree/$Format:%h$
%          $Format:%ci$

if ( ~isfield(ACEfit_Par,'nPerm') || isempty(ACEfit_Par.nPerm) )
    ACEfit_Par.nPerm = 1000;
end
if ( ~isfield(ACEfit_Par,'nBoot') || isempty(ACEfit_Par.nBoot) )
    ACEfit_Par.nBoot = 1000;
end
if ( ~isfield(ACEfit_Par,'nParallel') || isempty(ACEfit_Par.nParallel) )
    ACEfit_Par.nParallel = 1;
end
nPerm     = ACEfit_Par.nPerm;
nBoot     = ACEfit_Par.nBoot;
nParallel = ACEfit_Par.nParallel;

ACEfit_Par.RunID = 1:nParallel;

nMZF  = ACEfit_Par.nMZF;
nDZF  = ACEfit_Par.nDZF;
nSibF = ACEfit_Par.nSibF;
nSG   = ACEfit_Par.nSG;

% Permutation preparation
if nPerm>0
    
    if rem(nPerm,nParallel)
        error('''nPerm'' cannot be evenly divided by ''nParallel''!!')
    end
    
    ACEfit_Par.nPermPerRun = ones(nParallel,1)*nPerm/nParallel;
    
    if ( nMZF==0 || nDZF==0 )
        error('Cannot make permutation inference as there are no MZ or DZ twins!')
    end
    ACEfit_Par.Perm_index  = CreatePerm(nMZF,nDZF,nPerm);
    
end

% Bootstrap preparation
if nBoot>0
    
    if rem(nBoot,nParallel)
        error('''nBoot'' cannot be evenly divided by ''nParallel''!!')
    end
    
    ACEfit_Par.nBootPerRun = ones(nParallel,1)*nBoot/nParallel;

    ACEfit_Par.Boot_index  = CreateBoot(nMZF,nDZF,nSibF+nSG,nBoot);
    
end

% Save the main input file in the result folder
save(fullfile(ACEfit_Par.ResDir,'ACEfit_Par'),'ACEfit_Par');

return
