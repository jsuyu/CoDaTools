function results = mrlr_scores(X)
% mrlr_scores transforms composition to Log-Ratio scores centered by the midrange.
%
% midrange is the geometric mean of the maximum and the minimum components.
%
% Syntax:
%   results = mrlr_scores(X)
%
% Input:
%   X - Table or matrix (n × D) or  or vector (1 × D)
%       All values must be strictly positive (> 0)
%
% Output:
%   results - Table or structure with fields:
%               .scores       - mrlr scores (n × D) matrix or table
%               .scores_names - Names of mrlr scores (format: 'mrlr_i')
%
% Examples:
%   % Multiple observations
%   X = [0.3, 0.2, 0.3, 0.2;
%        0.25, 0.25, 0.25, 0.25];
%   results = mrlr(X);
% results.scores_names = {'mrlr_1', 'mrlr_2', 'mrlr_3', 'mrlr_4'}
%
% Author: Jordi
% Date: 2025
% Updated: 2025

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

    % Check that all values are positive
    [is_comp, row_check] = check_comp(X);
    if ~is_comp
        invalid_rows = find(~row_check);
        error( ['✗ Some rows contain zeros or negative values\n' ...
             '  Invalid rows: %d / %d\n' ...
             '  Row indices: %s'], ...
            numel(invalid_rows), n, mat2str(invalid_rows'));
    end
    
    % Apply log transformation
    X_log = log(X);
    
    % Center by subtracting row means
    row_midranges = (max(X_log,[], 2) + min(X_log,[], 2))/2;
    X_mrlr = X_log - row_midranges;

    %% Create coordinate names: alr_ij format
    scores_names = cell(1, D);
    for i = 1:D
        scores_names{i} = sprintf('mrlr_%d', i);
    end

    %% Package results
    if X_is_table
        % Return a table if input X was a table
        results = array2table(X_mrlr, 'VariableNames', scores_names);
    else
        % Return a struct if input X was an array
        results = struct( ...
            'scores_names', {scores_names}, ...
            'scores', X_mrlr ...
        );
    end

end
