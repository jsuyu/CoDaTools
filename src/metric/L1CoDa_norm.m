function Y = L1CoDa_norm(X)
% L1CoDa_norm computes induced L1 norm of a composition
%
% Calculates L1CoDa norm of compositional data
% The L1CoDa norm is the L1 Euclidean norm in mlr scores
%
% Syntax:
%   Y = L1CoDa_norm(X)
%
% Input:
%   X - Table or matrix (n × D) or  or vector (1 × D)
%       All values must be strictly positive (> 0)
%
% Output:
%   Y - norm vector
%
% Formula:
%   L1CoDa_norm(x) = ||mlr(x)||_1
%
% Example:
%   X = [0.25, 0.25, 0.25, 0.25;
%        0.40, 0.30, 0.20, 0.10;
%        0.10, 0.20, 0.30, 0.40];
%   Y = L1CoDa_norm(X);
%
% Author: Jordi
% Date: 2025
% Updated: 2025

    % Transform X to mlr scores
    X_mlr = mlr_scores(X);  % (n × D)
    
    if istable(X)
        Y = sum(abs(X_mlr),2);
        Y.Properties.VariableNames{1}='L1CoDa_norm';
    else
        Y =sum(abs(X_mlr.scores),2);  
    end
end