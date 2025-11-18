function Y = CoDa_distance(X1, X2, norm_type, option)
% CoDa_distance computes compositional distances between datasets
%
% Syntax:
%   D = CoDa_distance(X1)
%   D = CoDa_distance(X1, X2)
%   D = CoDa_distance(X1, X2, norm_type)
%   D = CoDa_distance(X1, X2, norm_type, option)
%
% Input:
%   X1        - Table or matrix (n1 × D) or vector (1 × D)
%               All values must be strictly positive (> 0)
%   X2        - (optional) Table or matrix (n2 × D) or vector
%               If empty or omitted, computes distances within X1
%   norm_type - (optional) Distance type:
%               'L1Coda', 'L1plr', 'L1clr', 'Ait', 'Linf'
%               Default: 'Ait'
%   option    - (optional) 'paired' or 'all-to-all'
%               Default: 'paired' if n1==n2, else 'all-to-all'
%               If X2=[], automatically uses 'all-to-all'
%
% Output:
%   Y - Distance matrix or vector
%       - If X2=[]: symmetric matrix (n1 × n1) with zeros on diagonal
%       - If option='paired': vector (n × 1)
%       - If option='all-to-all': matrix (n1 × n2)
%       - Table if X1 is table, numeric array otherwise
%
% Examples:
%   % Distances within X1
%   X1 = [0.25, 0.25, 0.25, 0.25;
%         0.40, 0.30, 0.20, 0.10];
%   D = CoDa_distance(X1);  % 2×2 symmetric matrix
%
%   % Paired distances
%   X2 = [0.50, 0.30, 0.15, 0.05;
%         0.30, 0.30, 0.30, 0.10];
%   D = CoDa_distance(X1, X2, 'Ait', 'paired');  % 2×1 vector
%
%   % With vectors
%   x1 = [0.25, 0.25, 0.25, 0.25];
%   x2 = [0.50, 0.30, 0.15, 0.05];
%   D = CoDa_distance(x1, x2, 'Ait', 'paired');  % scalar
%
% Author: Jordi
% Date: 2025
% Updated: November 2025

    %% Parse inputs
    if nargin < 2 || isempty(X2)
        X2 = [];
    end
    if nargin < 3 || isempty(norm_type)
        norm_type = 'Ait';
    end
    if nargin < 4
        option = '';
    end

    norm_type = lower(char(norm_type));
    option = lower(char(option));

    %% Validate norm type
    valid_norms = {'l1coda', 'l1plr', 'l1clr', 'ait', 'linf'};
    if ~ismember(norm_type, valid_norms)
        error('CoDa_distance:InvalidNorm', ...
            'norm_type must be one of: %s', strjoin(valid_norms, ', '));
    end
    
    %% Validate option
    valid_options = {'paired', 'all-to-all', ''};
    if ~isempty(option) && ~ismember(option, valid_options)
        error('CoDa_distance:InvalidOption', ...
            'option must be ''paired'' or ''all-to-all''');
    end

    %% Handle X1 - convert to array and store names
    
    % Handle column vector - transpose to row
    X1_is_table = istable(X1);
    if ~X1_is_table && isvector(X1) && size(X1, 1) > 1
        X1 = X1';
    end

    if X1_is_table
        X1_names = X1.Properties.RowNames;
        X1 = table2array(X1);
    else
        X1_names = [];
    end
    
    [n1, D1] = size(X1);

    if isempty(X1_names)
        X1_names = compose("Row%d", 1:n1);
    end
      
    % Check X1 is compositional
    [is_comp, row_check] = check_comp(X1, []);
    if ~is_comp
        invalid_rows = find(~row_check);
        error('CoDa_distance:X1NotCompositional', ...
            ['✗ X1 contains non-positive values\n' ...
             '  Invalid rows: %d / %d\n' ...
             '  Row indices: %s'], ...
            numel(invalid_rows), n1, mat2str(invalid_rows'));
    end

    %% Handle X2
    compute_within_X1 = isempty(X2);
    
    if compute_within_X1
        % Compute distances within X1
        X2 = X1;
        X2_names = X1_names;
        X2_is_table = X1_is_table;
        n2 = n1;
        D2 = D1;
    else
        X2_is_table = istable(X2);
        
        % Only transpose if it's an array (not table) AND column vector
        if ~X2_is_table && isvector(X2) && size(X2, 1) > 1
            X2 = X2';
        end
        
        if X2_is_table
            X2_names = X2.Properties.RowNames;
            X2 = table2array(X2);
        else
            X2_names = [];
        end
        
        [n2, D2] = size(X2);
        
        % Always ensure we have names
        if isempty(X2_names)
            X2_names = compose("Row%d", 1:n2);
        end
        
        % Check X2 is compositional
        [is_comp, row_check] = check_comp(X2, []);
        if ~is_comp
            invalid_rows = find(~row_check);
            error('CoDa_distance:X2NotCompositional', ...
                ['✗ X2 contains non-positive values\n' ...
                 '  Invalid rows: %d / %d\n' ...
                 '  Row indices: %s'], ...
                numel(invalid_rows), n2, mat2str(invalid_rows'));
        end
    end

    %% Check dimensions match
    if D1 ~= D2
        error('CoDa_distance:DimMismatch', ...
            'X1 and X2 must have same number of columns. Got %d and %d', D1, D2);
    end

    %% Determine option if not specified
    %% Determine option if not specified
    if isempty(option)
        if n1 == n2 && ~compute_within_X1
            option = 'paired';
        else
            option = 'all-to-all';
        end
    end
    
    %% Select norm function
    switch norm_type
        case 'l1coda'
            norm_func = @(X) L1CoDa_norm(X);
        case 'l1plr'
            norm_func = @(X) L1plr_norm(X);
        case 'l1clr'
            norm_func = @(X) L1clr_norm(X);
        case 'ait'
            norm_func = @(X) Ait_norm(X);
        case 'linf'
            norm_func = @(X) LinfCoDa_norm(X);
    end

    %% Compute distances
    switch option
        case 'paired'
            if n1 ~= n2
                error('CoDa_distance:SizeMismatch', ...
                    'For paired option, X1 and X2 must have same number of rows. Got %d and %d', n1, n2);
            end
            
            % Paired distances: norm(x1_i / x2_i)
            D = zeros(n1, 1);
            for i = 1:n1
                ratio = X1(i, :) ./ X2(i, :);
                D(i) = norm_func(ratio);
            end

        case 'all-to-all'
            % All pairwise distances
            D = zeros(n1, n2);
            for i = 1:n1
                for j = 1:n2
                    ratio = X1(i, :) ./ X2(j, :);
                    D(i, j) = norm_func(ratio);
                end
            end
    end

    %% Return table if X1 or X2 was a table
    if X1_is_table || X2_is_table
        if strcmp(option, 'paired')
            D = array2table(D, 'VariableNames', {'Distance'}, 'RowNames', X1_names);
        else
            D = array2table(D, 'VariableNames', X2_names, 'RowNames', X1_names);
        end
    end
  
    Y = D;
end
