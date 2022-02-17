function [Type,CIFTIext] = FileType(fn)
%
% Returns one of "NIFTI" or "CIFTI", depending on the format of the input
% file.
%
% For now, just uses file extension, with .img or .nii to detect NIFTI;
% all else is called CIFTI.  This is *very* rudimentary!  Only really can
% tell between these two formats, and, e.g., would very easily be tricked
% by a Analyze image.
%
% For convenience, if a numeric variable is passed, it is classified as "MATRIX"
%_______________________________________________________________________
% Version: http://github.com/NISOx-BDI/APACE/tree/$Format:%h$
%          $Format:%ci$

Type     = 'UNKOWN';
CIFTIext = '';

if ischar(fn)
    [pth,nam,ext,num] = spm_fileparts(lower(fn));
    if ~isempty(num)
        warning(['Volume/frame specified with comman notation ignored: ' fn])
    end
    if strcmp(ext,'.img')
        Type = 'NIFTI';
    elseif strcmp(ext,'.nii')
        ext2 = lower(spm_str_manip(fn(1:end-4),'e'));
        if ( strcmp(ext2,'dscalar') || strcmp(ext2,'dlabel') || strcmp(ext2,'dtseries') || strcmp(ext2,'pconn') || strcmp(ext2,'pscalar') )
            Type     = 'CIFTI';
            CIFTIext = ext2;
        else
            Type = 'NIFTI';
        end
    elseif strcmp(ext,'.gii')
        Type = 'GIFTI';
    end
elseif isnumeric(fn)
    Type = 'MATRIX';
    
end

return
