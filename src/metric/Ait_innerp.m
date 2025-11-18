function Y = Ait_innerp(X1, X2, option)
% Ait_innerp computes Aitchison inner product between compositions
% The Aitchison inner product is computed in clr subspace
%
% Syntax:
%   Y = Ait_innerp(X1)
%   Y = Ait_innerp(X1, X2)
%   Y = Ait_innerp(X1, X2, option)
%
% Input:
%   X1     - Table or matrix (n1 × D) or vector (1 × D)
%            All values must be strictly positive (> 0)
%   X2     - (optional) Table or matrix (n2 × D) or vector
%            If empty or omitted, computes inner products within X1
%   option - (optional) String: 'paired' or 'all-to-all'
%            - 'paired': compute X1(i,:) · X2(i,:) for each i
%            - 'all-to-all': compute X1(i,:) · X2(j,:) for all i,j
%            Default: 'paired' if n1==n2, else 'all-to-all'
%            If X2=[], automatically uses 'all-to-all'
%
% Output:
%   Y - Inner product result:
%       - If X2=[]: symmetric matrix (n1 × n1) with norms² on diagonal
%       - If option='paired': vector (n × 1)
%       - If option='all-to-all': matrix (n1 × n2)
%       - Table if X1 or X2 is table, numeric array otherwise
%
% Formula:
%   <x1, x2>_A = sum(clr(x1) .* clr(x2))
%   where clr(x) = log(x) - mean(log(x))
%
% Examples:
%   % Inner products within X1
%   X1 = [0.25, 0.25, 0.25, 0.25;
%         0.40, 0.30, 0.20, 0.10];
%   Y = Ait_innerp(X1);  % 2×2 symmetric matrix
%
%   % Paired inner products
%   X2 = [0.50, 0.30, 0.15, 0.05;
%         0.30, 0.30, 0.30, 0.10];
%   Y = Ait_innerp(X1, X2, 'paired');  % 2×1 vector
%
%   % All-to-all inner products
%   Y = Ait_innerp(X1, X2, 'all-to-all');  % 2×2 matrix
%
% Author: Jordi
% Date: 2025
% Updated: November 2025

    %% Parse inputs
    if nargin < 2 || isempty(X2)
        X2 = [];
    end
    if nargin < 3
        option = '';
    end
    
    option = lower(char(option));
    
    %% Validate option
    valid_options = {'paired', 'all-to-all', ''};
    if ~isempty(option) && ~ismember(option, valid_options)
        error('Ait_innerp:InvalidOption', ...
            'option must be ''paired'' or ''all-to-all''');
    end

    %% Handle X1
    X1_is_table = istable(X1);
    
    % Only transpose if it's an array (not table) AND column vector
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
    
    % Always ensure we have names
    if isempty(X1_names)
        X1_names = compose("Row%d", 1:n1);
    end
    
    % Check X1 is compositional
    [is_comp, row_check] = check_comp(X1, []);
    if ~is_comp
        invalid_rows = find(~row_check);
        error('Ait_innerp:X1NotCompositional', ...
            ['✗ X1 contains non-positive values\n' ...
             '  Invalid rows: %d / %d\n' ...
             '  Row indices: %s'], ...
            numel(invalid_rows), n1, mat2str(invalid_rows'));
    end
    
    % Transform X1 to clr
    X1_log = log(X1);
    X1_clr = X1_log - mean(X1_log, 2);
    
    %% Handle X2
    compute_within_X1 = isempty(X2);
    
    if compute_within_X1
        % Compute inner products within X1
        X2_clr = X1_clr;
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
        [is_comp, row_check] = CoDaTools.check.check_comp(X2, []);
        if ~is_comp
            invalid_rows = find(~row_check);
            error('Ait_innerp:X2NotCompositional', ...
                ['✗ X2 contains non-positive values\n' ...
                 '  Invalid rows: %d / %d\n' ...
                 '  Row indices: %s'], ...
                numel(invalid_rows), n2, mat2str(invalid_rows'));
        end
        
        % Transform X2 to clr
        X2_log = log(X2);
        X2_clr = X2_log - mean(X2_log, 2);
    end
    
    %% Check dimensions match
    if D1 ~= D2
        error('Ait_innerp:DimensionMismatch', ...
            'X1 and X2 must have same number of columns. Got %d and %d', D1, D2);
    end
    
    %% Determine option if not specified
    if isempty(option)
        if n1 == n2 && ~compute_within_X1
            option = 'paired';
        else
            option = 'all-to-all';
        end
    end
    
    %% Compute inner products
    switch option
        case 'paired'
            if n1 ~= n2
                error('Ait_innerp:SizeMismatch', ...
                    'For paired option, X1 and X2 must have same number of rows. Got %d and %d', n1, n2);
            end
            
            % Paired inner products
            Y = sum(X1_clr .* X2_clr, 2);  % (n × 1)
            
            if X1_is_table || X2_is_table
                Y = array2table(Y, 'VariableNames', {'Ait_innerp'}, 'RowNames', X1_names);
            end
            
        case 'all-to-all'
            % All pairwise inner products
            Y = X1_clr * X2_clr';  % (n1 × n2)
            
            if X1_is_table || X2_is_table
                Y = array2table(Y, 'VariableNames', X2_names, 'RowNames', X1_names);
            end
    end
end