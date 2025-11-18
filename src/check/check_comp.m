function [is_composition, row_check] = check_comp(X, display)
% chec_comp checks if a matrix contains strictly positive values
%
% A valid composition has all values strictly positive (> 0)
%
% Syntax:
%   [is_composition, row_check] = check_comp(X)
%   [is_composition, row_check] = check_comp(X, display)
%
% Input:
%   X        - Table or numeric matrix (n × D) or vector (1 × D)
%              Each row should contain only positive values.
%   display    - (optional) Logical flag to control printed output:
%                   true  → always display results
%                   false → suppress display
%                Default: true only if called without output arguments
%
% Outputs:
%   is_composition - TRUE if ALL rows have strictly positive values
%   row_check - Logical vector (n × 1) indicating which rows are valid
%
% Examples:
%   X = [0.5, 0.3, 0.2; 0.4, 0.4, 0.2];
%   [is_comp, row_check] = check_comp(X, true);  % Display info
%
%   X = [0.5, 0, 0.2];  % Contains zero
%   [is_comp, ~] = check_comp(X);
%
%   X = [0.5, 0.3, -0.1];  % Contains negative
%   [is_comp, ~] = check_comp(X);  % FALSE
%
% Author: Jordi
% Date: 2025
% Updated:2025

    % Validate display 
    if nargin < 2 || isempty(display)
        display = (nargout == 0);  % auto-display only if no outputs
    end
    if ~islogical(display) || ~isscalar(display)
        error('check_comp:InvalidDisplayFlag', ...
              'Display flag must be a logical scalar (true/false).');
    end

    % Validate input X
    if istable(X)
        X = table2array(X);
    end
    if ~isnumeric(X)
        error('check_comp:InvalidInput', 'Input must be a numeric matrix or vector');
    end
    
    % If column vector, convert to row
    if ~istable(X) && isvector(X) && size(X, 1) > 1
        X = X';
    end
    
    % Dimensions
    [n, D] = size(X);
    
    % CHECK: All elements must be strictly positive (> 0)
    row_check = all(X > 0, 2);  % TRUE for each row if all elements > 0
    
    % GLOBAL RESULT: TRUE if ALL rows are valid
    is_composition = all(row_check);
    
    % DISPLAY INFORMATION (only if no output is assigned)
    if display
        if is_composition
            fprintf('✓ All %d rows contain strictly positive values\n', n);
        else
            fprintf('✗ Some rows contain zeros or negative values\n');
            invalid_rows = find(~row_check);
            fprintf('  Invalid rows: %d / %d\n', numel(invalid_rows), n);
            fprintf('  Row indices: %s\n', mat2str(invalid_rows'));
        end
    end
end