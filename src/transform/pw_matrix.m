function H = pw_matrix(D)
% pw_matrix constructs pairwise difference matrix
%
% Creates matrix H for pairwise log-ratios from D components
%
% Syntax:
%   H = pw_matrix(D)
%
% Input:
%   D - Number of compositional parts/components
%
% Output:
%   H - Pairwise matrix (D x n_pairs), where n_pairs = D*(D-1)/2
%       Each column represents one pairwise contrast: part_i - part_j
%
% Example:
%   H = build_pairwise_matrix(4);
%   H is 4x6, columns correspond to all pairwise contrasts
%
% Author: Jordi
% Date: 2025
% Updated: 2025

    % Validate input
    if ~isscalar(D) || D < 2 || D ~= floor(D)
        error('build_pairwise_matrix:InvalidD', 'D must be an integer >= 2');
    end
    
    % Calculate number of pairs
    n_pairs = D * (D - 1) / 2;
    
    % Initialize matrix
    H = zeros(D, n_pairs);
    
    % Fill matrix with pairwise contrasts
    k = 1;
    for i = 1:D-1
        for j = (i+1):D
            H(i, k) = 1;   % +1 for component i
            H(j, k) = -1;  % -1 for component j
            k = k + 1;
        end
    end
end