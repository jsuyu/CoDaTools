function Y = Ait_norm(X)
% Ait_norm computes Aitchison (L2) norm of a composition
% The Aitchison norm is the Euclidean norm in clr subspace
%
% Syntax:
%   Y = Ait_norm(X)
%
% Input:
%   X - Table or matrix (n × D) or  or vector (1 × D)
%       All values must be strictly positive (> 0)
%
% Output:
%   Y - Norm vector (n x 1)
%
% Formula:
%   Ait_norm(x) = ||clr(x)||_2 = sqrt(sum(clr(x)^2))
%   where clr(x) = log(x) - mean(log(x))
%
% Examples:
%   % With array
%   X = [0.25, 0.25, 0.25, 0.25;
%        0.40, 0.30, 0.20, 0.10];
%   Y = Ait_norm(X);
%
%   % With table
%   data = readtable('compositions.xlsx');
%   X = data(:, {'Part1', 'Part2', 'Part3', 'Part4'});
%   Y = Ait_norm(X);
%
% Author: Jordi
% Date: 2025
% Updated:2025

    % Transform X to clr
    X_clr = clr_scores(X);  % (n × D)
    
    if istable(X)
        Y = sqrt(sum(X_clr.^2,2));
        Y.Properties.VariableNames{1}='Ait_norm';
    else
        Y =sqrt(sum(X_clr.scores.^2,2)); 
    end
end
