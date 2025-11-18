function Y = L1clr_norm(X)
% Calculates L1clr_norm of compositional data
% The L1clr_norm is the L1 Euclidean norm in clr scores
%
% Syntax:
%   Y = L1clr_norm(X)
%
% Input:
%   X - Table or matrix (n × D) or  or vector (1 × D)
%       All values must be strictly positive (> 0)
%
% Output:
%   Y - norm vector (n x 1)
%
% Formula:
%   L1clr_norm(x) = ||clr(x)||_1
%
% Example:
%   X = [0.25, 0.25, 0.25, 0.25;
%        0.40, 0.30, 0.20, 0.10;
%        0.10, 0.20, 0.30, 0.40];
%   Y = L1clr_norm(X);
%
% Author: Jordi
% Date: 2025
% Updated: 2025
    
    % Transform X to clr
    X_clr = clr_scores(X);  % (n × D)
    
    if istable(X)
        Y = sum(abs(X_clr),2);
        Y.Properties.VariableNames{1}='L1clr_norm';
    else
        Y =sum(abs(X_clr.scores),2); 
    end
end