function results = alr_coord(X, denom_idx)
% alr_coord transforms composition to Additive Log-Ratio coordinates
%
% Applies the alr (Additive Log-Ratio) transformation using a specified
% component as denominator. Returns structure with coordinates, and coordinate names.
%
% Syntax:
%   results = alr_coord(X, denom_idx)
%
% Input:
%   X - Table or matrix (n × D) or  or vector (1 × D)
%       All values must be strictly positive (> 0)
%   denom_idx - Index of denominator component (1 to D)
%               The component used as reference in log-ratios
%
% Output:
%   results - Table or structure with fields:
%               .coord       - alr coordinates (n × D-1) matrix or table
%               .coord_names - Names of alr coordinates (format: 'alr_ij')
%
% Formula:
%   alr_ij = log(x_i / x_j) where j is denom_idx
%
% Examples:
%   % With array (denominator index 4)
%   X = [0.3, 0.2, 0.3, 0.2];
%   results = alr_coord(X, 4);
%   % results.coord_names = {'alr_14', 'alr_24', 'alr_34'}
%
%   % With table (denominator at index 4)
%   X = array2table([0.3, 0.2, 0.3, 0.2], ...
%       'VariableNames', {'sleep', 'Sedentary', 'Light', 'MVPA'});
%   results = alr_coord(X, 4);
%   % results.Properties.VariableNames = {'alr_14', 'alr_24', 'alr_34'}
%
% Author: Jordi
% Date: 2025
% Updated: 2025

    %% Validate denominator index
    if nargin < 2
        error('alr_coord:MissingDenominator', 'Denominator index is required');
    end
    
    if ~isscalar(denom_idx) || denom_idx ~= floor(denom_idx)
        error('alr_coord:InvalidDenominator', 'Denominator index must be an integer');
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

    [n, D] = size(X);

    if denom_idx < 1 || denom_idx > D
        error('alr_coord:DenominatorOutOfRange', ...
              'Denominator index must be between 1 and %d (number of parts), got %d', ...
              D, denom_idx);
    end

    %% Check composition
    [is_comp, row_check] = check_comp(X);
    
    if ~is_comp
        invalid_rows = find(~row_check);
        error( ['✗ Some rows contain zeros or negative values\n' ...
             '  Invalid rows: %d / %d\n' ...
             '  Row indices: %s'], ...
            numel(invalid_rows), n, mat2str(invalid_rows'));
    end
    
    %% Create ALR contrast matrix Psi (D × D-1)
    I = eye(D-1);
    Psi = [I(1:(denom_idx-1), :); -ones(1, D-1); I(denom_idx:end, :)];
       
    %% Compute ALR coordinates
    X_alr = log(X)*Psi;

    %% Create coordinate names: alr_ij format
    coord_names = cell(1, D);
    for i = 1:D
        coord_names{i} = sprintf('alr_%d%d', i, denom_idx);
    end
    coord_names(denom_idx)=[];

    %% Package results
    if X_is_table
        % Return a table if input X was a table
        results = array2table(X_alr, 'VariableNames', coord_names);
    else
        % Return a struct if input X was an array
        results = struct( ...
            'coord_names', {coord_names}, ...
            'coord', X_alr ...
        );
    end
end