function [v] = coeffTable2kronCoeffs(C)
%kronCoeffs2coeffTable For a 2D function, convert from Kronecker
% coefficients to the matrix/table of polynomial coefficients
% representation.
%
%   Input: v - cell array of polynomial coefficients 
%        
%   Output: C - matrix/table of polynomial coefficients

n=2; 
degree = length(C)/2; 
v = cell(1,degree); 

% Arrange energy function coefficients from C into v
for k=2:degree
    v{k} = zeros(n^k,1);
    % Construct matrix idx where each row is the multi-index for one element of X
    idx = sort(tt_ind2sub(ones(1, k) * n, (1:n ^ k)'),2);
    [~,IA,~] = unique(idx,'rows');
    
    for i = IA.'
        alpha = sum(idx(i,:) == 1);
        beta = sum(idx(i,:) == 2);

        v{k}(i) = C(alpha+1, beta+1);
    end

    v{k} = kronMonomialSymmetrize(v{k},n,k);

end


end