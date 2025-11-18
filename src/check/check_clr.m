function [is_clr, row_check] = check_clr(X_clr, tolerance, display)
% check_clr checks if data are valid CLR (Centered Log-Ratio) scores
%
% clr scores must sum to zero for each observation (row)
%
% Syntax:
%   [is_clr, row_check] = check_clr(clrX)
%   [is_clr, row_check] = check_clr(clrX, tolerance)
%   [is_clr, row_check] = check_clr(clrX, tolerance, display)
%
% Input:
%   clrX      - Numeric matrix or table (n × D) or vector (1 × D)
%               Each row should sum to zero
%   tolerance - (optional) Maximum deviation from zero (default: 1e-12)
%   display   - (optional) Logical flag to print messages
%                   true  → always display
%                   false → suppress display
%               Default: true only if called without output arguments
%
% Outputs:
%   is_clr - TRUE if ALL rows sum to zero (within tolerance)
%   row_check - Logical vector (n × 1) indicating valid rows
%
% Examples:
%   % Valid clr scores
%   X = [1, -0.5, -0.5; 0.3, 0.3, -0.6];
%   [is_clr, row_check] = check_clr(X);  % is_clr = 1
%
%   % Invalid (doesn't sum to zero)
%   X = [1, 2, 3];
%   [is_clr, ~] = check_clr(X);  % is_clr = 0
%
% Author: Jordi
% Date: 2025
% Updated: 2025

    % Validate tolerance and default
    if nargin < 2 || isempty(tolerance)
    tolerance = 1e-12;
    end
    
    if nargin < 3 || isempty(display)
        display = (nargout == 0);  % auto-display if no outputs
    end
    
    if ~isnumeric(tolerance) || ~isscalar(tolerance) || tolerance < 0
        error('check_clr:InvalidTolerance', 'Tolerance must be a non-negative scalar.');
    end
    
    if ~islogical(display) || ~isscalar(display)
        error('check_clr:InvalidDisplay', 'Display must be a logical scalar (true/false).');
    end
    
    % Validate input
    if istable(X_clr)
    X_clr = table2array(X_clr);
    end
    
    if ~isnumeric(X_clr)
        error('check_clr:InvalidInput', 'Input must be numeric');
    end
    
    % Handle vector: convert column to row
    if ~istable(X_clr) && isvector(X_clr) && size(X_clr, 1) > 1
        X_clr = X_clr';
    end
    
    [n, D] = size(X_clr);
    
    % Calculate row sums
    row_sums = sum(X_clr, 2);
    
    % Check: all row sums ≈ 0 (within tolerance)
    row_check = abs(row_sums) < tolerance;
    
    % Global result
    is_clr = all(row_check);
    
    % Display if no output assigned
    if display
        if is_clr
            fprintf('✓ All %d rows sum to zero (tolerance: %.2e)\n', n, tolerance);
        else
            fprintf('✗ Some rows do not sum to zero (tolerance: %.2e)\n', tolerance);
            invalid_rows = find(~row_check);
            fprintf(' Invalid rows: %d / %d\n', numel(invalid_rows), n);
            fprintf(' Row indices: %s\n', mat2str(invalid_rows'));
        end
    end
end