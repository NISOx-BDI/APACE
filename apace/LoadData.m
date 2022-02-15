function [Y,YM] = LoadData(Vs)
%
%  Loads the data (actually, mostly already loaded)
%

switch Vs.FileType
    
    case 'NIFTI'
        
        if Vs.OneFile
            % Note that .dat is memory-mapped; it only gets loaded into RAM here
            Y = reshape(double(Vs.Im.dat(:)),[prod(Vs.Dim) Vs.n]);
        else
            if isempty(Vs.ContSel)
                Y = zeros([prod(Vs.Dim) Vs.n]);
                for i=1:Vs.n
                    Y(:,i) = reshape(double(Vs.Im(i).dat(:)),[prod(Vs.Dim) 1]);
                end
            else
                % Now we can assume 4D images, only 1 volume of which is to be selected
                Y = zeros([prod(Vs.Dim) Vs.n]);
                for i=1:Vs.n
                    Y(:,i) = reshape(double(Vs.Im(i).dat(:,:,:,Vs.ContSel(1))),[prod(Vs.Dim) 1]);
                end
            end
        end
        if isempty(Vs.ImMsk)
            YM = [];
        else
            YM = reshape(Vs.ImMsk.dat(:),[prod(Vs.Dim) 1]);
        end
        
    case 'CIFTI'
        
        if Vs.OneFile
            Y = reshape(double(Vs.Im.cdata),[prod(Vs.Dim) Vs.n]);
        else
            if isempty(Vs.ContSel)
                Y = zeros([prod(Vs.Dim) Vs.n]);
                for i=1:Vs.n
                    Y(:,i) = reshape(double(Vs.Im{i}.cdata),[prod(Vs.Dim) 1]);
                end
            else
                % cdata can be any dimension, always final dimension will be selected
                tmp    = size(Vs.Im{1}.cdata);
                FinDim = tmp(end);
                Y = zeros([prod(Vs.Dim) Vs.n]);
                for i=1:Vs.n
                    tmp    = reshape(double(Vs.Im{i}.cdata),[prod(Vs.Dim) FinDim]);
                    Y(:,i) = tmp(:,Vs.ContSel(1));
                end
            end
        end
        if isempty(Vs.ImMsk)
            YM = [];
        else
            if isnumeric(Vs.ImMsk)
                YM = reshape(Vs.ImMsk,[prod(Vs.Dim) 1]);
            else % must be CIFTI
                YM = reshape(Vs.ImMsk.cdata,[prod(Vs.Dim) 1]);
            end
        end
        
    case 'MATRIX'
        
        Y = reshape(Vs.Im,[prod(Vs.Dim) Vs.n]);
        
        if isempty(Vs.ImMsk)
            YM = [];
        else
            YM = reshape(Vs.ImMsk,[prod(Vs.Dim) 1]);
        end
        
end

return
