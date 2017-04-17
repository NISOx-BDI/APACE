function d = Dims(X)
% Generic support function
% For a single input, returns a vector of length Ndim(X), the dimensions of X.
% For NIFTI, input can be an array, and then a matrix with 7 columns is
% returned, corresponding to the maximal dimension information in the
% NIFTI header.
% For CIFTI, input can be an cell array, and then a matrix with 7 columns is
% returned with the cooresponding dimensional information.
%_______________________________________________________________________
% Version: http://github.com/nicholst/APACE/tree/$Format:%h$
%          $Format:%ci$

if isa(X,'nifti')
    if length(X)==1
        d = double(X.hdr.dim(2:(Ndim(X)+1)));
    else
        d = zeros(length(X),7);
        for i = 1:length(X)
            d(i,:) = double(X.hdr.dim(2:8));
        end
    end
elseif ( isa(X,'gifti') && isfield(X,'cdata') )
    d = size(X.cdata);
    if ( length(d)==2 && d(2)==1 )
        d = d(1);
    end
elseif ( iscell(X) && isa(X{1},'gifti') && isfield(X{1},'cdata') )
    d = ones(length(X),7);
    for i = 1:length(X)
        tmp = size(X{i}.cdata);
        d(i,1:length(tmp)) = tmp;
    end
elseif isnumeric(X)
    d = size(X);
    if ( length(d)==2 && d(2)==1 )
        d = d(1);
    end
else
    error('Unknown image type')
end

return
