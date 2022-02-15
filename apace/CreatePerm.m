function Perm_index = CreatePerm(nG1,nG2,nPerm)
%
% Create the permutation relablings
%
% nG1    - number of subjects or subject pairs in group 1
% nG2    - number of subjects or subject pairs in group 2
% nPerm  - number of permutations
%
if nPerm>nchoosek(nG1+nG2,nG1)-1
    error('Too many permutations requested! Please assign %d or a smaller number to ''nPerm''!',nchoosek(nG1+nG2,nG1)-1)
end
% if nPerm>mynchoosek(nG1+nG2,nG1)
%     error('Too many permutations requested! Please set %d or a smaller number for ''nPerm''!',mynchoosek(nG1+nG2,nG1))
% end

Perm_index = zeros(nPerm,nG1+nG2);

% Derive seed from time
rng('shuffle');

for i=1:nPerm
    while(1)
        rand_label = randperm(nG1+nG2);
        G1_label   = rand_label(1:nG1);
        G2_label   = rand_label((nG1+1):end);
        Perm_label = [sort(G1_label) sort(G2_label)];
        if ( Perm_label(nG1)>nG1 && ~any(all(Perm_index(1:(i-1),:)==kron(ones(i-1,1),Perm_label),2)) )
            Perm_index(i,:) = Perm_label;
            break;
        end
    end
end

return

function n = mynchoosek(n1,n2)

n = exp(gammaln(n1+n2+1)-gammaln(n1+1)-gammaln(n2+1));

return
