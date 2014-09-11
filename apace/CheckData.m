function Vs = CheckData(ACEfit_Par)
%
% Vs  Structure describing key aspects of input image data
%
% Vs.n          - The number of subjects
% Vs.OneFile    - Flags use of a single file containing all subjects' data
% Vs.Dim        - The dimensions of the data for a single subject; note, in the
%                 case Vs.OneFile==1 and the data being supplied as a single k-D
%                 (e.g. 4D) image, Vs.Dim describes the (k-1)-dimensions
%                 (e.g. 3D) for each subject (not the original specified dataset).
% Vs.ClustInf   - Flags 3D images per subject, meaning cluster inf. possible
% Vs.ContSel    - "Contrast Selection", usable only in two cases:
%                  a. A list of 4D nifti files specified, and ContSel
%                     selects out the volume/frame to use.
%                  b. A list of cifit files specified, and ContSel
%                     selects out of the last dimension.
% Vs.ExampleImg - Information on a template image; used for writing
%                 images in the same format as read.
% Vs.Im         - Information on image data (possibly just a handle or actual
%                 data; depends on format; actual data reading takes place in
%                 LoadData.m).
% Vs.ImMks      - Information on mask data.
%

%
% Fill-in defaults
%
if ~isfield(ACEfit_Par,'Pmask');   ACEfit_Par.Pmask   = ''; end
if ~isfield(ACEfit_Par,'Dsnmtx');  ACEfit_Par.Dsnmtx  = ''; end
if ~isfield(ACEfit_Par,'Nlz');     ACEfit_Par.Nlz     = 1;  end
if ~isfield(ACEfit_Par,'ContSel'); ACEfit_Par.ContSel = []; end
if ~isfield(ACEfit_Par,'NoImg');   ACEfit_Par.NoImg   = 0;  end


P_nm  = ACEfit_Par.P_nm;
Pmask = ACEfit_Par.Pmask;

%%%%  Path must be configured to include wb_command (or make this an
%%%%  absolute path); ciftiopen and ciftisave must also be on Matlab path
WBC = 'wb_command';

% Check output results directory
[stat,msg] = mkdir(ACEfit_Par.ResDir);  % Need 2nd output arg to suppress warnings
if stat~=1
    error('Cannot create results directory')
end
if mkdir(ACEfit_Par.ResDir,'__TestDir')~=1
    error('Results directory is not writable')
else
    rmdir(fullfile(ACEfit_Par.ResDir,'__TestDir'))
end

% Number of all subjects
n    = size(ACEfit_Par.kin,1);
Vs.n = n;  % redundant, but useful

%
% Check family data and load the data as well (because you have to load
% it to check it).
%

Vs.OneFile  = 0;  % Flags a single file containing all subjects' data
Vs.ClustInf = 0;  % Flags 3D images per subject, meaning that
% spatial/cluster infernece is possible
Vs.ContSel  = []; % Possibly select a single image/row

if isstruct(P_nm)
    error('Direct reading of structures not supported')
    % Maybe could add something where V from spm_vol could be
    % specified... but not yet.
end


if ( iscell(P_nm) || ischar(P_nm) )
    
    %
    % Data on disk
    %
    
    % Extract file paths
    if iscell(P_nm)
        % P_nm is cell string of file paths
        if length(P_nm)==1
            Vs.OneFile = 1;
        end
        P = P_nm;
    elseif ischar(P_nm)
        % Load a text file of filepaths
        fpP = fopen(P_nm);
        if fpP<0
            error('Cannot open text filelist of images');
        end
        P = textscan(fpP,'%[^\n]');
        fclose(fpP);
        P = P{1};
    end
    
    if ( ~Vs.OneFile && length(P)~=n )
        error('Number of image filepaths specified doesn''t match kinship table');
    end
    if ( Vs.OneFile && ~isempty(ACEfit_Par.ContSel) )
        error('Cannot do contrast selection when using a single image file');
    end
    
    if strcmp(FileType(P{1}),'NIFTI')
        
        if isempty(Pmask)
            NM = [];
        else
            NM = nifti(Pmask);
        end
        
        if Vs.OneFile
            N = nifti(P{1});
            if N.hdr.dim(Ndim(N)+1)~=n
                % UNKWNON/TO-TEST: Can I trust Nifti's dimension (N.hdr.dim(1)==2) to *always* flag 2D files?
                disp(N.hdr.dim);
                error('Single input file''s last dimension doesn''t match kinship table (NIFTI)');
            end
            if Ndim(N)==4
                % This "OneFile" NIFTI is 4D; spatial inference is possible
                Vs.ClustInf = 1;
            end
            if ~isempty(NM)
                if Ndim(NM)~=Ndim(N)-1
                    error('Mask is wrong dimension; must be same dimension as each subject''s data')
                end
                tmp = Dims(N);
                if ~all(Dims(NM)==tmp(1:end-1))
                    error('Mask and each subject''s data have different dimensions')
                end
            end
            tmp    = Dims(N);
            Vs.Dim = tmp(1:end-1);
            
        else
            
            N = nifti(P);
            if any(any(diff(Dims(N))))
                error('Subjects have inconsistent image dimensions (NIFTI)')
            end
            if Ndim(N(1))==3
                % Each NIFTI is 3D; spatial inference is possible
                Vs.ClustInf = 1;
            end
            if ~isempty(ACEfit_Par.ContSel)
                if Ndim(N(1))~=4
                    error('Cannot do contrast selection unless NIFTI files are 4D');
                else
                    Vs.ContSel = ACEfit_Par.ContSel;
                end
                tmp    = Dims(N(1));
                Vs.Dim = tmp(1:end-1);
            else
                Vs.Dim = Dims(N(1));
            end
            
            if ~isempty(NM)
                if ( Ndim(NM)~=length(Vs.Dim) || ~all(Dims(NM)==Vs.Dim) )
                    error('Mask doesn''t match dimension of each subject''s data (NIFTI)')
                end
            end
        end
        
        Vs.ExampleImg = N(1);   % Used by clone_vol
        Vs.FileType   = 'NIFTI';
        Vs.Im         = N;
        Vs.ImMsk      = NM;
        
    elseif strcmp(FileType(P{1}),'CIFTI')
        
        if CheckCIFTIfun(WBC)==0
            error('No CIFTI support')
        end
        
        if Vs.OneFile
            C   = ciftiopen(P{1},WBC);
            tmp = Dims(C);
            if tmp(end)~=n
                disp(size(C));
                error('Single input file''s last dimension doesn''t match kinship table (CIFTI)');
            end
            Vs.Dim              = tmp(1:end-1);
            Vs.ExampleImg       = C;
            Vs.ExampleImg.cdata = [];  % Save some RAM
        else
            C = cell(1,0);
            for i=1:n
                C{i} = ciftiopen(P{i},WBC);
            end
            if any(any(diff(Dims(C))))
                error('Subjects have inconsistent image dimensions (CIFTI)')
            end
            
            if ~isempty(ACEfit_Par.ContSel)
                Vs.ContSel = ACEfit_Par.ContSel;
                tmp        = Dims(C{1});
                Vs.Dim     = tmp(1:end-1);
            else
                Vs.Dim     = Dims(C{1});
            end
            Vs.ExampleImg  = C{1};
        end
        
        if isempty(Pmask)
            CM = [];
        else
            if strcmp(FileType(Pmask),'MATRIX')
                CM = Pmask;
            elseif strcmp(FileType(Pmask),'CIFTI')
                CM = ciftiopen(Pmask,WBC);
            else
                error('CIFTI images can only be masked with CIFTI or matrix data');
            end
            if ( Ndim(CM)~=length(Vs.Dim) || ~all(Dims(CM)==Vs.Dim) )
                error('Mask doesn''t match dimension of each subject''s data (CIFTI)')
            end
        end
        
        [~,CIFTIext] = FileType(P{1});
        Vs.CIFTIext  = CIFTIext;
        Vs.FileType  = 'CIFTI';
        Vs.Im        = C;
        Vs.ImMsk     = CM;
    else
        error('Unrecognized file type')
    end
    
else
    
    %
    % Data already in Matlab workspace
    %
    
    if ~isempty(ACEfit_Par.ContSel)
        error('Cannot do contrast selection on in-memory matrix--do it yourself!');
    end
    
    if isempty(Pmask)
        XM = [];
    else
        XM = Pmask;
    end
    
    X = P_nm;
    if size(X,Ndim(X))~=n
        error('Data last dimension doesn''t match kinship table (MATRIX)');
    end
    if Ndim(X)==4
        % Could be wrong, but seems reasonable to assume the first 3 dims are space
        Vs.ClustInf = 1;
    end
    if ~isempty(XM)
        if Ndim(XM)~=Ndim(X)-1
            error('Mask is wrong dimension; must be same dimension as each subject''s data')
        end
        tmp = Dims(X);
        if ~all(Dims(XM)==tmp(1:end-1))
            error('Mask and each subject''s data have different dimensions')
        end
    end
    
    Vs.FileType = 'MATRIX';
    tmp         = Dims(X);
    Vs.Dim      = tmp(1:end-1);
    Vs.Im       = X;
    Vs.ImMsk    = XM;
    
end

return

