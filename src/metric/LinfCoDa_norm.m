function Y = LinfCoDa_norm(X)
% LinfCoDa_norm computes induced Linf norm of a composition
%
% Calculates LinfCoDa norm of compositional data
% The LinfCoDa norm is the Linf Euclidean norm in mrlr scores
%
% Syntax:
%   Y = LinfCoDa_norm(X)
%
% Input:
%   X - Table or matrix (n × D) or  or vector (1 × D)
%       All values must be strictly positive (> 0)
%
% Output:
%   Y - norm vector
%
% Formula:
%   LinfCoDa_norm(x) = ||mrlr(x)||_inf
%
% Example:
%   X = [0.25, 0.25, 0.25, 0.25;
%        0.40, 0.30, 0.20, 0.10;
%        0.10, 0.20, 0.30, 0.40];
%   Y = LinfCoDa_norm(X);
%
% Author: Jordi
% Date: 2025
% Updated: 2025
    
     % Transform X to mrlr
    X_mrlr = mrlr_scores(X);  % (n × D)
    
    if istable(X)
        Y = max(abs(X_mrlr),[],2);
        Y.Properties.VariableNames{1}='L1CoDa_norm';
    else
        Y = max(abs(X_mrlr.scores),[],2);
    end
end