function Psi = pv_basis(D, order)
% pv_basis creates orthonormal pivot basis
%
% Constructs an orthonormal basis using sequential binary partitions
% where each component is compared against the remaining ones
%
% Syntax:
%   Psi = pivot_basis(D)
%   Psi = pivot_basis(D, order)
%
% Input:
%   D - Number of compositional parts
%   order - (optional) Order of components (default: [1,2,...,D])
%           Permutation vector specifying the pivot sequence
%
% Output:
%   Psi - Orthonormal pivot basis (D × D-1)
%         Column i: component order(i) vs remaining components
%
% Example:
%   % Default order [1,2,3,4]
%   Psi = pivot_basis(4);
%   % Column 1: {1} vs {2,3,4}
%   % Column 2: {2} vs {3,4}
%   % Column 3: {3} vs {4}
%
%   % Custom order [2,1,3,4]
%   Psi = pivot_basis(4, [2,1,3,4]);
%   % Column 1: {2} vs {1,3,4}
%   % Column 2: {1} vs {3,4}
%   % Column 3: {3} vs {4}
%
% Author: Jordi
% Date: 2025

    % Validate D
    if ~isscalar(D) || D < 2 || D ~= floor(D)
        error('pivot_basis:InvalidD', 'D must be an integer >= 2');
    end
    
    % Default order
    if nargin < 2
        order = 1:D;
    end
    
    % Validate order
     if ~isvector(order) || length(order) ~= D || ~all(sort(order) == 1:D)
        error('pv_basis:InvalidOrder', 'Order must be a permutation of 1:%d', D);
     end
    
    % Initialize basis
    Psi = zeros(D, D-1);
    
    % Build each column
    for col = 1:D-1
        pivot = order(col);
        neg_idx = order(col+1:end);  % remaining components after pivot

        r = 1;                  % positive group size
        s = length(neg_idx);    % negative group size

        coef_pos = sqrt(s / (r * (r + s)));
        coef_neg = -sqrt(r / (s * (r + s)));

        Psi(pivot, col) = coef_pos;
        Psi(neg_idx, col) = coef_neg;
    end
end