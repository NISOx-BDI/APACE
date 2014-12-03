function ACEfit_Par = Reorder(ACEfit_Par)
%
% Reorder the input data
%
% InfMx - Kinship information matrix; a CSV file with a header row and then
%         one row per subject. There must be 4 columns: an subject ID
%         (ignored), a mother ID, a father ID, and zygosity.
%
% Possible types of relationships in the kinship matrix:
%
% -1    - Self
%  0    - Unrelated
%  1    - MZ twins
%  2    - DZ twins
%  3    - Siblings
%

%
% Read mother ID, father ID and zygosity
%
fid  = fopen(ACEfit_Par.InfMx);
if fid==-1
    error('Cannot open kinship information CSV file')
end
Afid = textscan(fid,'%s','Delimiter', ',');
fclose(fid);

n    = length(Afid{1})/4 - 1;
Bfid = reshape(Afid{1},4,n+1)';
Bfid = Bfid(2:end,:);
Bfid = cellfun(@(x) x(x~=' ' & x~='"'),Bfid,'un',0);

FamIDs = cell2mat(cellfun(@str2num,strcat(Bfid(:,2),Bfid(:,3)),'un',0));

Zyg = zeros(n,1);
Zyg(strcmp('MZ',     Bfid(:,4))) = 1;
Zyg(strcmp('NotMZ',  Bfid(:,4))) = 2;
Zyg(strcmp('NotTwin',Bfid(:,4))) = 3;


%
% Arrange subjects from the same family together
%
[sFamIDs,IDs] = sort(FamIDs);
Zyg           = Zyg(IDs);


%
% Rearrange families in order of MZ families, DZ families, non-twin 
% families, and singletons
%
uFamIDs = unique(sFamIDs,'stable');
FamZyg  = inf(size(sFamIDs));
i0      = 0;

% Initialize the number of MZ, DZ, and non-twin families
nMZF  = 0;
nDZF  = 0;
nSibF = 0;

for i=1:length(uFamIDs)
    
    nF = sum(sFamIDs==uFamIDs(i));
    
    if nF>1
        iZyg = Zyg(sFamIDs==uFamIDs(i));
        
        if sum(iZyg)==2*min(iZyg)+(nF-2)*3
            FiZyg = min(iZyg);
            
            if FiZyg==1
                nMZF  = nMZF+1;
            elseif FiZyg==2
                nDZF  = nDZF+1;
            else
                nSibF = nSibF+1;
            end
            
        else
            FiZyg = 3;
            nSibF = nSibF+1;
        end
        
        FamZyg(i0+[1:nF]) = FiZyg;
    end
    
    i0 = i0+nF;
    
end

[~,FamZygIDs] = sort(FamZyg);
IDs           = IDs(FamZygIDs);
sFamIDs       = sFamIDs(FamZygIDs);
Zyg           = Zyg(FamZygIDs);


%
% For each twin family, twins are put first
%
uFamIDs = unique(sFamIDs,'stable');
nFam    = zeros(size(uFamIDs));
Tmp     = zeros(n,1);
i0      = 0;

for i=1:length(uFamIDs)
    
    nF             = sum(sFamIDs==uFamIDs(i));
    nFam(i)        = nF;
    
    [~,iF]         = sort(Zyg(i0+[1:nF]));
    TmpF           = i0+[1:nF];
    Tmp(i0+[1:nF]) = TmpF(iF);
    
    i0             = i0+nF;
    
end

IDs = IDs(Tmp);
Zyg = Zyg(Tmp);

nSG = sum(nFam==1);

%
% Generate kinship matrix after reordering
%
[MZMx,DZMx,SibMx] = deal(zeros(n));

for fam=1:length(uFamIDs)
    
    Ifam  = double(sFamIDs==uFamIDs(fam));
    SibMx = SibMx + Ifam*Ifam';
    
    Ifam  = double((sFamIDs==uFamIDs(fam)) & (Zyg==1));
    MZMx  = MZMx + Ifam*Ifam';
    
    Ifam  = double((sFamIDs==uFamIDs(fam)) & (Zyg==2));
    DZMx  = DZMx + Ifam*Ifam';
    
end

kin = SibMx*3 - MZMx*2 - DZMx;
kin = kin - diag(diag(kin)) - eye(n);

%
% Update input information
%
ACEfit_Par.nMZF  = nMZF;   % Number of MZ families
ACEfit_Par.nDZF  = nDZF;   % Number of DZ families
ACEfit_Par.nSibF = nSibF;  % Number of non-twin families
ACEfit_Par.nSG   = nSG;    % Number of singletons

ACEfit_Par.nFam  = nFam;   % Number of subjects for each family
ACEfit_Par.kin   = kin;    % Kinship matrix
ACEfit_Par.ID    = IDs;    % Reordered ID's

fprintf(['*** Data summary\n',...
   	'  %4d Subjects in total\n',...
	'  %4d MZ families\n',...
	'  %4d DZ families\n',...
	'  %4d Non-twin families\n',...
	'  %4d Singletons\n',...
	'  %4d Smallest family size\n',...
	'  %4d Largest family size\n\n'],...
	sum(nFam),nMZF,nDZF,nSibF,nSG,min(nFam),max(nFam));


return

