function WriteData(X,Vs,fname,fpath)
%
%  Write data, based on template image
%
%_______________________________________________________________________
% Version: http://github.com/nicholst/APACE/tree/$Format:%h$
%          $Format:%ci$

if length(Vs.Dim)==1
    Dim = [Vs.Dim 1];
else
    Dim = Vs.Dim;
end

switch Vs.FileType
    
    case 'NIFTI'
        
        N     = Vs.ExampleImg;
        N.dat = file_array(fullfile(fpath,[fname '.nii']),Dim,...
            [spm_type('float32') spm_platform('bigend')]);
        create(N); % File is now mapped; setting .dat writes to disk
        
        % This seems idiotic, but it seems necessary to use the mapped file class
        switch length(Dim)
            case 2
                N.dat(:,:) = reshape(X,Dim);
            case 3
                N.dat(:,:,:) = reshape(X,Dim);
            case 4
                N.dat(:,:,:,:) = reshape(X,Dim);
            case 5
                N.dat(:,:,:,:,:) = reshape(X,Dim);
            case 6
                N.dat(:,:,:,:,:,:) = reshape(X,Dim);
            case 7
                N.dat(:,:,:,:,:,:,:) = reshape(X,Dim);
        end
        
    case 'CIFTI'
        
        C       = Vs.ExampleImg;
        C.cdata = reshape(X,Dim);
        ciftisave(C,fullfile(fpath,[fname '.' Vs.CIFTIext '.nii']),'wb_command');
        
    case 'MATRIX'
        
        % Save a variable named fname to a mat file named fname
        tmp = struct(fname,reshape(X,Dim));
        save(fullfile(fpath,fname),'-struct','tmp');
        
end

return
