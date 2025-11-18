function Y = LpCoDa_norm(X, p)
% LpCoDa_norm computes induced LpCoDa norm of a composition
%
% Computes the Lp norm of compositional data, defined as:
%   LpCoDa_norm(x) = min_lambda ||clr(x) + lambda * 1_D||_p
%
% Syntax:
%   Y = LpCoDa_norm(X, p)
%
% Input:
%   X - Table or matrix (n × D) or  or vector (1 × D)
%       All values must be strictly positive (> 0)
%   p - Norm order (p >= 1)
%
% Output:
%   Y - Norm vector (n × 1)
%
% Example:
%   X = [0.25, 0.25, 0.25, 0.25;
%        0.40, 0.30, 0.20, 0.10;
%        0.10, 0.20, 0.30, 0.40];
%   Y = LpCoDa_norm(X, 2);
%
% Author: Jordi
% Date: 2025
% Updated: 2025

    if nargin < 2
        p = 2;  % Default to L2CoDa (Aitchison norm)
    end
    if ~isscalar(p) || p < 1
        error('LpCoDa_norm:InvalidP', 'p must be a scalar >= 1');
    end
    
    X_is_table=false;
    %% Parse input: table or array
    if istable(X)
        X_is_table=true;
        X = table2array(X);
    end

    % Transform X to CLR
    X_clr = clr_scores(X).scores;  % (n × D)
    
    [n, D] = size(X_clr);
    
    % Initialize output
    Y = zeros(n,1);
    
    % Define ones vector
    onesD = ones(1, D);
    
    % Loop over compositions
    for i = 1:n
        z = X_clr(i, :);
        % Objective function in lambda
        f = @(lambda) sum(abs(z + lambda*onesD).^p);
        % Minimize over lambda (real number)
        L = max(abs(z));  % safe bound for lambda
        [~, f_opt] = fminbnd(f, -L, L);  % adjust range if needed
        % Compute minimized value
        Y(i) = (f_opt)^(1/p);
    end

    if X_is_table
        Y = array2table(Y);
        Y.Properties.VariableNames{1}='LpCoDa_norm';
    end
end
