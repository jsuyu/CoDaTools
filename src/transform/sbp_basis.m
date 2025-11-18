function Psi = sbp_basis(SBP, tolerance)
% sbp_basis computes orthonormal olr basis from SBP matrix
%
% Validates SBP matrix and computes the orthonormal basis
%
% Syntax:
%   Psi = sbp_basis(SBP)
%   Psi = sbp_basis(SBP, tolerance)
%
% Input:
%   SBP - table or matrix (D x (D-1))
%         Must contain only -1, 0, +1
%   tolerance - (optional) Orthogonality tolerance (default: 1e-12)
%
% Output:
%   Psi - Orthonormal basis (D x (D-1))
%
% Checks performed (throws error if fails):
%   1. Dimensions
%   2. Values: only -1, 0, +1
%   3. Partitions: each row has +1 and -1
%   4. Orthogonality: Psi'*Psi = I
%
% Example:
%   SBP = [ 1,  1,  0;
%           1, -1,  0;
%          -1,  0,  1;
%          -1,  0, -1];
%   Psi = sbp_basis(SBP);
%
% Author: Jordi
% Date: 2025
% Updated: 2025

    % Validate tolerance
    if nargin < 2
        tolerance = 1e-12;
    end
    
    if ~isnumeric(tolerance) || ~isscalar(tolerance) || tolerance < 0
        error('sbp_basis:InvalidTolerance', 'Tolerance must be a non-negative scalar');
    end

    % Validate input
    if istable(SBP)
        SBP = table2array(SBP);
    end

    % Validate SBP using check_sbp
    check_sbp(SBP, tolerance);

    % Get dimensions
    [D , K] = size(SBP);
    
    % Compute orthonormal basis Psi
    Psi = zeros(D, K);
    
    for i = 1:K
        % Count positive and negative components
        r = sum(SBP(:, i) == 1);
        s = sum(SBP(:, i) == -1);
        
        % Compute coefficients
        coef_pos = sqrt(s / (r * (r + s)));
        coef_neg = -sqrt(r / (s * (r + s)));
        
        % Fill basis vector
        Psi(SBP(:, i) ==  1, i) = coef_pos;
        Psi(SBP(:, i) == -1, i) = coef_neg;
    end
end