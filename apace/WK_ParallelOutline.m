%
% Preparation
%

% % Uncomment this to hide figures (for the rest of your Matlab session!).
% % This is useful if you don't want to see the large number of figure
% % windows that get created.
% set(0, 'DefaultFigureVisible', 'off');


%%% 1) Create a structure array specifying details of data
ACEfit_Par.Model   = 'ACE';             % Choose a model (AE or ACE) for 
                                        % data fitting.

ACEfit_Par.P_nm    = 'APACE_imgs.txt';  % A list of images or other data
                                        % specification; see FileFormats.txt
                                        
ACEfit_Par.InfMx   = 'KinInf.csv';      % Kinship information matrix of 4
                                        % columns with headers:
                                        % 'SubjectID', 'MotherID',
                                        % 'FatherID', and 'Zygosity'.
                                        % MotherID and FatherID must be
                                        % numeric, and Zygosity is 'MZ',
                                        % 'NotMZ' or 'NotTwin'. The order
                                        % of subjects must match that of
                                        % image paths or subjects in
                                        % ACEfit_Par.P_nm.
ACEfit_Par.ResDir  = '/my/path/ResDir';

%%% The rest are optional; omit to use default values

ACEfit_Par.Pmask   = 'APACE_mask.nii';  % Brain mask image (default: whole volume)

ACEfit_Par.Dsnmtx  = '';                % Design matrix (default: all-ones vector)

ACEfit_Par.Nlz     = 1;                 % Inverse Gaussian normalisation options 
                                        % (applied to each phenotype voxel/element)
                                        %  0 - None
                                        %  1 - Gaussian normalisation *before* 
                                        %      forming residuals
                                        %  2 - Gaussian normalisation *after* 
                                        %      forming residuals; note, though, 
                                        %      that after this Gaussian 
                                        %      normalisation the residuals
                                        %      will be formed again to ensure
                                        %      orthognality of residuals to
                                        %      the design matrix.  

ACEfit_Par.AggNlz  = 0;                 % Aggregate heritability normalisation options
                                        % (applied to each phenotype voxel/element)
                                        %  0 - Default, meaning that only
                                        %      determined by Nlz option.
                                        %      This always entails
                                        %      de-meaning and removal of any effects
                                        %      in Dsnmtx, but variance at each
                                        %      voxel/element unchanged.  
                                        %  1 - Same as 0, but with variance
                                        %      normalisation; this is full
                                        %      'studentization'.
                                        %  2 - *Undo* mean centering of
                                        %      residual formation. (If Dsnmtx
                                        %      is empty or default value of
                                        %      an all-ones vector, this is
                                        %      equivalent to using the raw
                                        %      input data. If Dsnmtx
                                        %      contains nuisance regressors
                                        %      the residuals will be formed
                                        %      and the mean added back in. 
                                        %  3 - Same as 2, but with
                                        %      variance normalisation;
                                        %      note that for each
                                        %      voxel/element, mean is
                                        %      first added back and then
                                        %      data are divided by stdev.

ACEfit_Par.ContSel = [];                % Select a single contrast (a volume
                                        % in a 4D Nifti file supplied for
                                        % each subject, or at the last
                                        % dimension in a cifti image; NOT
                                        % compatible with a single file
                                        % containing all subjects' data.)
                                        
ACEfit_Par.NoImg   = 0;                 % If 1, suppress image-wise inference,
                                        % and only compute summaries.
                                        
                                        
%%% 2) Updata 'ACEfit_Par' with the input data information
ACEfit_Par = PrepData(ACEfit_Par);

%%% 3) Run the original data once
ACEfit_Par.alpha_CFT = [];              % Cluster-forming threshold (default: 0.05)

ACEfit_Par = ACEfit(ACEfit_Par);
%%% 4) Add permutation and bootstrapping information, and save "ACEfit_Par.mat"
ACEfit_Par.nPerm = 1000;                % Number of permutations
ACEfit_Par.nBoot = 1000;                % Number of bootstrap replicates
nParallel        = [];                  % Number of parallel runs (default: 1, without parallelization)

PrepParallel(ACEfit_Par,nParallel);


%
% Permutation inference for computing FWE-corrected p-values
%

if ACEfit_Par.nPerm>0
    
    %%% 1)
    %%%%% The following code can be scripted or parallelized as you wish.  
    %%%%% Please refer to "README_APACE_intro.pdf" for the example snippets
    %%%%% for parallelization.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load(fullfile(ACEfit_Par.ResDir,'ACEfit_Par.mat'));
    RunID = 1;
    ACEfit_Perm_Parallel(ACEfit_Par,RunID);
    % Inside ResDir, it will create result sets, ACEfit_Parallel_XXXX.mat,
    % one for each RunID.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% 2) This will merge together all the results from the specified RunID's
    load(fullfile(ACEfit_Par.ResDir,'ACEfit_Par.mat'));
    ACEfit_Perm_Parallel_Results(ACEfit_Par);
    
    %%% 3)
    FWEalpha = 0.05;
    FDRalpha = 0.05;
    load(fullfile(ACEfit_Par.ResDir,'ACEfit_Par.mat'));
    ACEfit_Results(ACEfit_Par,FWEalpha,FDRalpha);
    
end


%
% Bootstrapping inference for constructing CIs
%

if ACEfit_Par.nBoot>0
    
    %%% 1)
    %%%%% The following code can be scripted or parallelized as you wish.
    %%%%% Please refer to "README_APACE_intro.pdf" for the example snippets
    %%%%% for parallelization.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load(fullfile(ACEfit_Par.ResDir,'ACEfit_Par.mat'));
    RunID = 1;
    ACEfit_Boot_Parallel(ACEfit_Par,RunID);
    % Inside ResDir, it will create result sets, BootCI_Parallel_XXXX.mat,
    % one for each RunID.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% 2) This will merge together all the results from the specified RunID's
    load(fullfile(ACEfit_Par.ResDir,'ACEfit_Par.mat'));
    ACEfit_Boot_Parallel_Results(ACEfit_Par);
    
    %%% 3)
    Balpha = 0.05;
    load(fullfile(ACEfit_Par.ResDir,'ACEfit_Par.mat'));
    Boot_CIs(ACEfit_Par,Balpha)
    
end


%
% Aggregate heritability (aka "Steve's method") for multiple phenotypes, 
% with P-values via permutation and CI's via boostrapping.
%
% Note that permuation (steps 1-3) and bootstrapping (steps 1-3) can be
% skipped by setting ACEfit_Par.nPerm=0 and ACEfit_Par.nBoot=0 separately; 
% the following code can be run immediately after "PrepParallel". 
%
load(fullfile(ACEfit_Par.ResDir,'ACEfit_Par.mat'));
Palpha = 0.05; % Significance threhold for permutations (plotting only)
Balpha = 0.05; % Confidence level, where CI's have level 100*(1-alpha)

AgHe_Method(ACEfit_Par,Palpha,Balpha)     
% % Once with no variance normalisation
% ACEfit_Par.AggNlz  = 0;   % de-meaning only
% AgHe_Method(ACEfit_Par,Palpha,Balpha,'_NoNorm')
% 
% % Now, again, but with Variance normalisation
% ACEfit_Par.AggNlz  = 1;   % de-meaning and scaling to have stdev of 1.0
% AgHe_Method(ACEfit_Par,Palpha,Balpha,'_Norm')


%
% Generate summary file 
%

APACEsummary(ACEfit_Par,'ResultSummary')


