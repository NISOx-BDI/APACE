function ACEfit_Par = ProcData(ACEfit_Par,Y,YM)
%
% Process data for further computation
%
%_______________________________________________________________________
% Version: http://github.com/nicholst/APACE/tree/$Format:%h$
%          $Format:%ci$

nElm = size(Y,1);

n      = ACEfit_Par.n;
Subset = ACEfit_Par.Subset;

Y = Y(:,Subset);

X = ACEfit_Par.Dsnmtx;
if isempty(X)
    X = ones(n,1);
elseif ischar(X)
    X = load(X);
    X = X(Subset,:);
else
    X = X(Subset,:);
end
if size(X,1)~=n
    error('Design matrix has wrong size! Expect %d rows, found %d. \n',n,size(X,1))
end

%
% Find the indices of MZ, DZ & Sib pairs
%
kin     = triu(ACEfit_Par.kin);
MZcode  = 1;
DZcode  = 2;
Sibcode = 3;
fk_MZ   = find(kin==MZcode);
fk_DZ   = find(kin==DZcode);
fk_Sib  = find(kin==Sibcode);
ACEfit_Par.I_MZ  = [rem(fk_MZ,n)  (fk_MZ  - rem(fk_MZ,n) )/n+1];
ACEfit_Par.I_DZ  = [rem(fk_DZ,n)  (fk_DZ  - rem(fk_DZ,n) )/n+1];
ACEfit_Par.I_Sib = [rem(fk_Sib,n) (fk_Sib - rem(fk_Sib,n))/n+1];

% Save original mean, if needed later
if ( ACEfit_Par.AggNlz==2 || ACEfit_Par.AggNlz==3 )
    ACEfit_Par.Ymean = mean(Y,2);
end

% Inverse Gaussian transform
if ACEfit_Par.Nlz>0
    EY    = norminv((1:n)/(n+1));
    EY    = EY/std(EY);
    pinvX = (X*pinv(X))';
    for j = 1:nElm
        if ACEfit_Par.Nlz==2
            Res    = Y(j,:) - Y(j,:)*pinvX;
            [~,oY] = sort(Res);
        else
            [~,oY] = sort(Y(j,:));
        end
        Y(j,oY) = EY*std(Y(j,:))+mean(Y(j,:));
    end
end

% Reorder design matrix
ACEfit_Par.X = X(ACEfit_Par.ID,:);
% Reorder the data
ACEfit_Par.Y = Y(:,ACEfit_Par.ID);

% In-mask data index after removing NaN's and null vectors
tmp = ones(nElm,1);
for j = 1:nElm
    tmp(j) = tmp(j) && ~any(isnan(Y(j,:))) && ~all(diff(Y(j,:))==0);
    if ~isempty(YM)
        tmp(j) = tmp(j) && YM(j);
    end
end
I_data = find(tmp);
if isempty(YM) && ~all(tmp==1)
    fprintf('Voxels/elements removed due to missingness or being constant (%d).\n',sum(~tmp));
elseif ~isempty(YM) && sum(~tmp)>sum(~YM)
    fprintf('Voxels/elements removed due to missingness or being constant (%d; plus %d due to mask).\n',sum(~tmp)-sum(~YM),sum(~YM));
end

if isempty(I_data)
    error('All in-mask data vectors are null vectors or contain NaN')
else
    ACEfit_Par.I_data = I_data';
end

% Check the selected model
if ( ~strcmpi(ACEfit_Par.Model,'ACE') && ~strcmpi(ACEfit_Par.Model,'AE') )
    error('Unknown model. The fitted model should be ''AE'' or ''ACE''! \n')
end

return
