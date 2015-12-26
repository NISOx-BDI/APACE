%
% Preparation
%

set(0, 'DefaultFigureVisible', 'off');

%%% 1) Create a structure array specifying details of data
ACEfit_Par.P_nm       = 'APACE_imgs.txt';
                                            % A list of images or other data
                                            % specification; see FileFormats.txt
                                            % The order of subjects must match 
                                            % the order of image paths or subjects 
                                            % in ACEfit_Par.InfMx.

ACEfit_Par.InfMx      = 'KinInf.csv';       % Kinship information matrix of 4
                                            % columns with headers
                                            % ('SubjectID','MotherID',
                                            % 'FatherID','Zygosity'). MotherID
                                            % and FatherID must be numeric,
                                            % and Zygosity is one of 'MZ',
                                            % 'NotMZ' and 'NotTwin'.
                                            
ACEfit_Par.ResDir     = '/my/path/ResDir';


%%% The rest are optional; omit to use default values

ACEfit_Par.Pmask      = 'APACE_mask';       % Brain mask image (default: whole volume)
    
ACEfit_Par.Dsnmtx     = '';                 % Design matrix (default: all-ones vector)

ACEfit_Par.Nlz        = 1;                  % Inverse Gaussian normalisation options
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
    
ACEfit_Par.ContSel    = [];                 % Select a single contrast (a volume
                                            % in a 4D Nifti file supplied for
                                            % each subject, or at the last
                                            % dimension in a cifti image; NOT
                                            % compatibile with a single file
                                            % containing all subjects' data.)
                                                                                                                            
   
%%% 2) Updata 'ACEfit_Par' with the input data information
ACEfit_Par       = PrepData(ACEfit_Par);

%%% 3) Run the original data once
ACEfit_Par       = ACEfit(ACEfit_Par);

%%% 4) Add permutation information, and save "ACEfit_Par.mat"
ACEfit_Par.nPerm = 1000;                    % Number of permutations (should be positive integer).
nParallel        = [];                      % Number of parallel runs (default: 1, without parallelization)

PrepParallel(ACEfit_Par,nParallel);

%
% Permutations
%

if ACEfit_Par.nPerm>0
    
    %%% 1) For each paired phenotypes, the following code can be scripted or parallelized.
    %%%%% Please refer to part (2A) of "README_APACE_intro.pdf" in univariate
    %%%%% APACE software for the example snippets for parallelization.
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

