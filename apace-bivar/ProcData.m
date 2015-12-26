function ACEfit_Par = ProcData(ACEfit_Par,Y,YM)
%
% Process data for further computation
%

[nElm,n] = size(Y);

X = ACEfit_Par.Dsnmtx;
if isempty(X)
    X = ones(n,1);
elseif ischar(X)
    X = load(X);
end
if size(X,1)~=n
    error('Design matrix has wrong size! Expect %d rows, found %d. \n',n,size(X,1))
end

% Inverse Gaussian transform
if ACEfit_Par.Nlz>0
    EY = norminv((1:n)/(n+1));
    EY = EY/std(EY);
    if ACEfit_Par.Nlz==2
        Res    = Y - Y*(X*pinv(X))';
        [~,oY] = sort(Res,2);
    else
        [~,oY] = sort(Y,2);
    end
    for j = 1:nElm
        Y(j,oY(j,:)) = EY*std(Y(j,:)) + mean(Y(j,:));
    end
end

% Reorder design matrix & data
ACEfit_Par.X = X(ACEfit_Par.ID,:);
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
% if isempty(YM)
%     I_data = ~all(diff(Y,1,2)==0,2) & ~any(isnan(Y),2);
% else
%     I_data = ~all(diff(Y,1,2)==0,2) & ~any(isnan(Y),2) & YM(:);
% end

if isempty(I_data)
    error('All in-mask data vectors are null vectors or contain NaN')
else
    ACEfit_Par.I_data = I_data';
end

return
