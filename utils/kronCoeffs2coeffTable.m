function [C] = kronCoeffs2coeffTable(v)
%kronCoeffs2coeffTable For a 2D function, convert from Kronecker
% coefficients to the matrix/table of polynomial coefficients
% representation.
%
%   Input: v - cell array of polynomial coefficients 
%        
%   Output: C - matrix/table of polynomial coefficients

n=2; 
degree = length(v);

% Arrange energy function coefficients from v into C matrix
Csize = degree+1;
C = zeros(Csize,Csize); % Just allocate plenty of space for C so that the code doesn't try to call and go outside of C


for k=2:degree
    % Construct matrix idx where each row is the multi-index for one element of X
    idx = tt_ind2sub(ones(1, k) * n, (1:n ^ k)');
    
    for i = 1:n^k
        alpha = sum(idx(i,:) == 1);
        beta = sum(idx(i,:) == 2);
        C(alpha+1, beta+1) = C(alpha+1, beta+1) + v{k}(i);
    end
end

end