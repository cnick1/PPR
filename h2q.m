function [q] = h2q(h)
%h2q Given a cell array of the coefficients for h(x), compute the
% coefficients of q(x) such that q(x) = h(x).' * h(x)
%   Detailed explanation goes here

% Create a vec function for readability
vec = @(X) X(:);

% Process inputs
n = size(h{1},2);
ell = length(h);

% Construct q(x) coefficients
for k = 2:2*ell
    q{k} = sparse(n^k,1);
    
    % TODO: use symmetry to cut in half
    for i = (k - ell):ell % would be 1:(k-1) but need to truncate only to h{} terms which exist
        j = k - i;
        q{k} = q{k} + vec(h{i}.' * h{j});
    end
end

end

