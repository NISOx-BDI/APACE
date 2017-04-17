%_______________________________________________________________________
% Version: http://github.com/nicholst/APACE/tree/$Format:%h$
%          $Format:%ci$

function OK = CheckCIFTIfun(WBC)
global CIFTI_OK
if ~isempty(CIFTI_OK)
    OK = CIFTI_OK;
else
    OK = 0;
    if exist('ciftiopen')==2 && exist('ciftisave')==2
        if WBC(1) == filesep
            % Absolute path
            if exist(WBC,'file')==2
                OK = 1;
            end
        else
            % OK... so this is totally Unix dependent
            if system(['which ' WBC])==0
                OK = 1;
            end
        end
    end
    CIFTI_OK = OK;
end

return
