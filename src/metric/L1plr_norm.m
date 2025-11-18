function Y = L1plr_norm(X)
% L1plr_norm computes the absolute log-ratio norm of compositions
%
% Computes:
%   norm(x) = (1 / (D - 1)) * sum_{i<j} |log(x_i / x_j)|
%
% Syntax:
%   Y = LmeanAbsLogRatio_norm(X)
%
% Input:
%   X - Table or matrix (n × D) or  or vector (1 × D)
%       All values must be strictly positive (> 0)
%
% Output:
%   Y - Norm vector (n × 1)
%
% Example:
%   X = [0.25, 0.25, 0.25, 0.25;
%        0.40, 0.30, 0.20, 0.10;
%        0.10, 0.20, 0.30, 0.40];
%   Y = L1plr_norm(X);
%
% Author: Jordi
% Date: 2025
% Updated: 2025

    % Handle vector: convert column to row
    if ~istable(X) && isvector(X) && size(X, 1) > 1
        X = X';
    end
    
    % Check composition
    [is_comp, row_check] = check_comp(X);
    if ~is_comp
        invalid_rows = find(~row_check);
        error( ['✗ Some rows contain zeros or negative values\n' ...
             '  Invalid rows: %d / %d\n' ...
             '  Row indices: %s'], ...
            numel(invalid_rows), n, mat2str(invalid_rows'));
    end
    
    [n , D] = size(X);

    if istable(X)
        X = table2array(X);
        X_pw = log(X)*pw_matrix(D); % Compute pairwise logratios (n × p), p = D*(D-1)/2
        Y = array2table(sum(abs(X_pw), 2)/(D - 1), 'VariableNames',{'L1plr_norm'});
    else
        X_pw = log(X)*pw_matrix(D); % Compute pairwise logratios (n × p), p = D*(D-1)/2
        Y =sum(abs(X_pw), 2)/(D - 1);
    end
end
