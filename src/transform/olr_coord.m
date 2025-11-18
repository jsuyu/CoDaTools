function results = olr_coord(X, Psi, tolerance)
% olr_coord transform composition to orthonormal Log-Ratio coordinates
%
% Applies olr transformation using an orthonormal basis
%
% Syntax:
%   results = olr(X, Psi)
%   results = olr(X, Psi, tolerance)
%
% Input:
%   X - Table or matrix (n × D) or  or vector (1 × D)
%       All values must be strictly positive (> 0)
%   Psi - Orthonormal basis (D × D-1)
%         Each column is a basis vector (e.g., from sbp_basis)
%   tolerance - (optional) Orthonormality tolerance (default: 1e-12)
%
% Output:
%   results - Table or structure with fields:
%               .coord       - olr scores (n × (D-1)) matrix or table
%               .coord_names - Names of mrlr scores (format: 'olr_i')
%
% Formula:
%   results = clr_scores(X) * Psi
%
% Example:
%   % Create orthonormal basis from SBP
%   SBP = [ 1,  1,  0;
%           1, -1,  0;
%          -1,  0,  1;
%          -1,  0, -1];
%   Psi = sbp_basis(SBP);
%   
%   % Transform composition
%   X = [0.35, 0.15, 0.20, 0.30];
%   results = olr(X, Psi);
%   results.coord_names = {'olr_1', 'olr_2', 'olr_3'}
%
% Author: Jordi
% Date: 2025
% Updated: 2025

     if nargin < 2 || isempty(Psi)
            error('olr_coord:MissingInput', ...
                  'At least two inputs (X and Psi) are required.');
     end

    % Default tolerance
    if nargin < 3
        tolerance = 1e-12;
    end
    
    if ~isnumeric(Psi)
        error('olr:InvalidBasis', 'Basis Psi must be numeric');
    end
    
    if ~isnumeric(tolerance) || ~isscalar(tolerance) || tolerance < 0
        error('olr:InvalidTolerance', 'Tolerance must be a non-negative scalar');
    end
    
    [D_Psi, K] = size(Psi);

    % Check that Psi is orthonormal
    max_error = max(abs(Psi' * Psi - eye(K)), [], 'all');
    if max_error >= tolerance
        error('olr:NotOrthonormal', ...
              'Basis Psi is not orthonormal (max error: %.2e, tolerance: %.2e)', ...
              max_error, tolerance);
    end

    X_is_table=false;
     %% Parse input: table or array
    if istable(X)
        X_is_table=true;
        X = table2array(X);
    end

    % Handle vector: convert column to row
    if ~X_is_table && isvector(X) && size(X, 1) > 1
        X = X';
    end

    [n , D] = size(X);
    
    % Validate dimensions
    if D ~= D_Psi
        error('olr:DimensionMismatch', ...
            'X has %d columns but Psi has %d rows', D, D_Psi);
    end
    
    if K ~= D - 1
        error('olr:InvalidBasisDimension', ...
            'Basis Psi must have D-1 columns');
    end

    % Transform to CLR
    X_clr = clr_scores(X).scores;  % (n × D)
    
    % Project onto orthonormal basis
    X_olr = X_clr * Psi;  % (n × D) * (D × D-1) = (n × D-1)

    %% Create coordinate names: alr_ij format
    coord_names = cell(1, D-1);
    for i = 1:(D-1)
        coord_names{i} = sprintf('olr_%d', i);
    end

    %% Package results
    if X_is_table
        % Return a table if input X was a table
        results = array2table(X_olr, 'VariableNames', coord_names);
    else
        % Return a struct if input X was an array
        results = struct( ...
            'coord_names', {coord_names}, ...
            'coord', X_olr ...
        );
    end
end