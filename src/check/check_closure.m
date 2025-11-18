function [is_closed, row_check, closure_value] = check_closure(X, tolerance, display)
% check_clousure checks if data have constant row sums (closure property)
%
% A closed composition has:
%   1. All values strictly positive (> 0)
%   2. All rows sum to the same constant
%
% Syntax:
%   [is_closed, row_check, closure_value] = check_closure(X)
%   [is_closed, row_check, closure_value] = check_closure(X, tolerance)
%   [is_closed, row_check, closure_value] = check_closure(X, tolerance, display)
%
% Input:
%   X - Table or matrix (n × D) or vector (1 × D)
%       Each row should be a composition with constant sum
%   tolerance - (optional) Maximum deviation from common sum (default: 1e-12)
%   display    - (optional) Logical flag to control printed output:
%                   true  → always display results
%                   false → suppress display
%                Default: true only if called without output arguments
%
% Outputs:
%   is_closed - TRUE if all rows are positive and sum to same value
%   row_check - Logical vector (n × 1) indicating valid rows
%   closure_value - Common sum value (or NaN if not constant)
%
% Examples:
%   % Closed to 1 (proportions)
%   X = [0.5, 0.3, 0.2; 0.4, 0.4, 0.2];
%   [is_closed, ~, k] = check_closure(X);  % is_closed=1, k=1
%
%   % Closed to 100 (percentages)
%   X = [50, 30, 20; 40, 40, 20];
%   [is_closed, ~, k] = check_closure(X);  % is_closed=1, k=100
%
%   % Force display even when outputs are captured 
%   [a,b,c] = check_closure(X, 1e-10, true);
%
% Author: Jordi
% Date: 2025
%Updated:2025

    % Validate tolerance and display
     if nargin < 2 || isempty(tolerance)
        tolerance = 1e-12;
    end

    if nargin < 3 || isempty(display)
        display = (nargout == 0);  % auto-display only if no output
    end

    if ~islogical(display) || ~isscalar(display)
        error('check_closure:InvalidDisplay', ...
              'Display must be a logical scalar (true/false).');
    end

    if ~isnumeric(tolerance) || ~isscalar(tolerance) || tolerance < 0
        error('check_closure:InvalidTolerance', ...
              'Tolerance must be a non-negative scalar.');
    end
    
    if istable(X)
        X = table2array(X);
    end

    % Handle vector: convert column to row
    if ~istable(X) && isvector(X) && size(X, 1) > 1
        X = X';
    end

    % Check composition (includes numeric validation and vector conversion)
    [is_composition, row_positive] = check_comp(X);
    [n,D]=size(X);

    % Stop if not compositional
    if ~is_composition
        is_closed = false;
        row_check = row_positive;
        closure_value = NaN;
        if display
            fprintf('✗ Some rows contain zeros or negative values\n');
            invalid_rows = find(~row_positive);
            fprintf('  Invalid rows: %d / %d\n', numel(invalid_rows), n);
            fprintf('  Row indices: %s\n', mat2str(invalid_rows'));
        end
        return;
    end
    
    % Calculate row sums
    row_sums = sum(X, 2);
    
    % Check if all sums are equal (within tolerance)
    reference_sum = mode(round(row_sums/tolerance)*tolerance);
    row_check = abs(row_sums - reference_sum) < tolerance;
    
    % Global result
    is_closed = all(row_check);
    
    % Determine closure value
    if is_closed
        closure_value = reference_sum;
    else
        closure_value = NaN;
    end
    
    % Display if no output assigned
    if display
        if is_closed
            fprintf('✓ All %d rows are closed to %.6f (tolerance: %.2e)\n', ...
                n, closure_value, tolerance);
        else
            invalid_rows = find(~row_check);
            fprintf('✗ Invalid rows: %d / %d (tolerance: %.2e)\n', ...
                numel(invalid_rows), n, tolerance);
            fprintf('  Expected sum: %.6f\n', reference_sum);
            fprintf('  Row indices: %s\n', mat2str(invalid_rows'));
        end
    end
end