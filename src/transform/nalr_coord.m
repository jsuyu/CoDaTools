function results = nalr_coord(X, index_sets)
% nalr_coord transform composition using nested (recursive) alr schema.
%
% Creates log-ratios recursively using nested sets of denominator components.
% Each set uses the geometric mean of its components as denominator.
%
% Syntax:
%   results = nalr_coord(X, index_sets)
%
% Input:
%   X - Table or matrix (n × D) or  or vector (1 × D)
%       All values must be strictly positive (> 0)
%   index_sets - Cell array of nested index sets {set1, set2, ..., setK}
%                Each set contains column indices to use as denominator
%                Rules:
%                - First set: at most D-1 components
%                - Last set: exactly 1 component
%                - Each set is strictly contained in the previous one
%                - At most D-1 sets
%
% Output:
%   results - Table or structure with fields:
%               .coord       - nalr scores (n × D) matrix or table
%               .coord_names - Names of mrlr scores (format: 'mrlr_i')
%       Number of columns from set i: length(set_i-1) - length(set_{i})
%
% Example:
%   % D=4 components
%   X = [0.2, 0.3, 0.25, 0.25];
%   
%   % Set 1: {3, 4} → log(X₁/geomean(X₃,X₄)), log(X₂/geomean(X₃,X₄))
%   % Set 2: {4}   → log(X₃/X₄)
%   index_sets = {[3, 4], [4]};
%   results = nalr_coord(X, index_sets);  % Returns 3 coordinates
%   results.coord_names = {'nalr_1/[ 3 4 ]', 'nalr_2/[ 3 4 ]', 'nalr_3/[ 4 ]'}
%
% Author: Jordi
% Date: 2025
% Updated: 2025

    if ~iscell(index_sets) || isempty(index_sets)
        error('nested_alr_coord:InvalidIndexSets', 'index_sets must be a non-empty cell array');
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

    % Check that all values are positive
    [is_comp, row_check] = check_comp(X);
    if ~is_comp
        invalid_rows = find(~row_check);
        error( ['✗ Some rows contain zeros or negative values\n' ...
             '  Invalid rows: %d / %d\n' ...
             '  Row indices: %s'], ...
            numel(invalid_rows), n, mat2str(invalid_rows'));
    end
    
    % Validate index_sets
    K = numel(index_sets);
    if K > D - 1
        error('nested_alr_coord:TooManySets', ...
            'Number of sets (%d) exceeds D-1 (%d)', K, D-1);
    end
    
    index_sets=[{1:D},index_sets];

    % Validate last set size
    last_set = index_sets{end};
    if numel(last_set) ~= 1
        error('nested_alr_coord:LastSetNotSingle', ...
            'Last set must contain exactly 1 component, got %d', numel(last_set));
    end

    % Validate each set
    for k = 2:K+1
        curr_set = index_sets{k};
        
        % Check it's a vector of integers
        if ~isvector(curr_set) || ~all(curr_set == floor(curr_set))
            error('nested_alr_coord:InvalidSet', ...
                'Set %d must be a vector of integers', k);
        end
        
        % Check indices are in valid range
        if any(curr_set < 1) || any(curr_set > D)
            error('nested_alr_coord:IndicesOutOfRange', ...
                'Set %d contains indices outside [1, %d]', k, D);
        end
        
        % Check for duplicates
        if length(unique(curr_set)) ~= length(curr_set)
            error('nested_alr_coord:DuplicateIndices', ...
                'Set %d contains duplicate indices', k);
        end
        
        % Validate nesting: each set must be contained in the previous one
        prev_set = index_sets{k-1};
        if ~all(ismember(curr_set, prev_set))
            error('nested_alr_coord:NotNested', ...
                'Set %d is not strictly contained in set %d', k, k-1);
        end
    end
    
    % Initialize output
    X_nalr = zeros(n, D - 1);
    coord_names = cell(1, D-1);
    col_idx = 1;  % Current column in X_nalr
    
    % Apply log transformation
    X_log = log(X);
    
    for k = 2:K+1
        prev_set=index_sets{k-1};
        curr_set = index_sets{k};
        
        % Components NOT in current set but in previous set
        numerator_idx = setdiff(prev_set, curr_set);
        
        % Geometric mean of denominator (components in current set)
        % In log-scale: mean of logs
        denom_log = mean(X_log(:, curr_set), 2);
        denom_str=['[ ' sprintf('%d ', curr_set) ']'];

        % Compute log-ratios for numerator components
        for j = numerator_idx
            X_nalr(:, col_idx) = X_log(:, j) - denom_log;
            coord_names{col_idx} = sprintf('nalr_%d/%s', j,denom_str);
            col_idx = col_idx + 1;
        end
    end

    %% Package results
    if X_is_table
        % Return a table if input X was a table
        results = array2table(X_nalr, 'VariableNames', coord_names);
    else
        % Return a struct if input X was an array
        results = struct( ...
            'coord_names', {coord_names}, ...
            'coord', X_nalr ...
        );
    end
end