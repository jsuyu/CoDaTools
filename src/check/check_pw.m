function [is_pw, X] = check_pw(pw, D, tolerance,display)
% check_pw checks if data are valid pairwise scores and recover composition
%
% Verifies that pw are valid pairwise log-ratios and recovers the 
% corresponding composition X
%
% Syntax:
%   [is_pw, D, X] = check_pw(pw)
%   [is_pw, D, X] = check_pw(pw, D)
%   [is_pw, D, X] = check_pw(pw, D, tolerance)
%   [is_pw, D, X] = check_pw(pw, D, tolerance, display)
%
% Input:
%   pw - Table or matrix (n × p) or vector of pairwise scores
%        where p = D*(D-1)/2
%   D - (optional) Number of compositional parts
%       If not provided, computed from size of pw
%   tolerance - (optional) Reconstruction error tolerance (default: 1e-12)
%   display   - (optional) Logical flag to print messages
%                   true  → always display
%                   false → suppress display
%               Default: true only if called without output arguments
%
% Outputs:
%   is_pw - TRUE if all rows are valid pairwise coordinates
%   D - Number of compositional parts (computed or provided)
%   X - Recovered composition matrix (n × D) with rows closed to 1
%       Returns NaN if not valid pairwise
%
% Examples:
%   % Auto-detect D
%   [is_pw, D, X] = check_pw(pw);  % D computed from size of pw
%
%   % Specify D explicitly
%   [is_pw, D, X] = check_pw(pw, 4);
%
% Author: Jordi
% Date: 2025
% Updated: 2025

    % Validate tolerance and display
    if nargin < 3 || isempty(tolerance)
        tolerance = 1e-12;
    end

    if nargin < 4 || isempty(display)
        display = (nargout == 0);  % auto-display if no outputs
    end
    
    if ~isnumeric(tolerance) || ~isscalar(tolerance) || tolerance < 0
        error('check_pw:InvalidTolerance', 'Tolerance must be a non-negative scalar.');
    end

    if ~islogical(display) || ~isscalar(display)
        error('check_pw:InvalidDisplay', 'Display must be a logical scalar (true/false).');
    end
    
    pw_is_table=false;
    % Validate input
    if istable(pw)
        pw_is_table=true;
        pw = table2array(pw);
    end

    if ~isnumeric(pw)
        error('check_pw:InvalidInput', 'Input pw must be numeric');
    end
    
    % Handle vector: convert column to row
    if ~pw_is_table && isvector(pw) && size(pw, 1) > 1
        pw = pw';
    end
    
    % Get dimensions
    [n, p] = size(pw);
    
    % STEP 1: Determine or validate D
    if nargin < 2 || isempty(D)
        % Compute D from p = D*(D-1)/2
        discriminant = 1 + 8*p;
        if discriminant < 0 || sqrt(discriminant) ~= floor(sqrt(discriminant))
            error('check_pw:InvalidDimensions', ...
                'Number of columns (%d) does not correspond to D*(D-1)/2 for any integer D', p);
        end
        D = (1 + sqrt(discriminant)) / 2;
        if D ~= floor(D)
            error('check_pw:InvalidDimensions', ...
                'Number of columns (%d) does not correspond to D*(D-1)/2 for any integer D', p);
        end
    else
        % Validate provided D
        if ~isscalar(D) || D < 2 || D ~= floor(D)
            error('check_pw:InvalidD', 'D must be an integer >= 2');
        end
        expected_p = D * (D - 1) / 2;
        if p ~= expected_p
            error('check_pw:DimensionMismatch', ...
                'For D=%d, expected %d columns but got %d', D, expected_p, p);
        end
    end
    
    % STEP 2: Build pairwise matrix H
    H = pw_matrix(D);
    
    % STEP 3: Compute pseudoinverse of H
    H_pinv = pinv(H);
    
    % STEP 4: Recover log-coordinates
    X_log = pw * H_pinv;  % (n x D)
    
    % STEP 5: Check reconstruction (verify pw is in range of H)
    pw_reconstructed = X_log * H;  % (n × p)
    errors = sqrt(sum((pw - pw_reconstructed).^2, 2));
    error_norm = max(errors);
    is_pw = error_norm < tolerance;
    
    % STEP 6: Recover composition from log-coordinates
    if is_pw
        % Apply exponential and closure
        X = exp(X_log);  % (n × D)
        row_sums = sum(X, 2);
        X = X ./ row_sums;  % Close to sum = 1
    else
        % Return NaN if not valid pairwise
        X = NaN(n, D);
    end
    
    coord_names = cell(1, D);
    for i = 1:D
        coord_names{i} = sprintf('c%d', i);
    end

    if pw_is_table
        X = array2table(X, 'VariableNames', coord_names);
    else
        X = struct( ...
            'comp', X, ...
            'variables', {coord_names} ...
        );
    end

    % Display if no output assigned
    if display
        if is_pw
            fprintf('✓ All %d rows are valid pairwise coordinates for D=%d\n', n, D);
            fprintf('  Max reconstruction error: %.2e (tolerance: %.2e)\n', ...
                error_norm, tolerance);
            fprintf('  Composition X recovered successfully\n');
            n_show = min(6, n);
            D_show = min(5, D);  % first 5 columns
            last_col = D;         % always include last column
            
            fprintf('  Showing first %d of %d rows and selected columns:\n', n_show, n);
            
            col_width = 10;  % width for each column
            
            if pw_is_table
                if D > 6
                    % Select first 5 and last column
                    T = X(:, [1:D_show, last_col]);
                    disp(T);
                    fprintf('    ... (columns %d to %d omitted) ...\n', D_show+1, D-1);
                else
                    disp(X(1:n_show, :));
                end
            else
                % Display variable names
                fprintf('    ');
                for j = 1:D_show
                    fprintf('%*s', col_width, X.variables{j});
                end
                if D > 6
                    fprintf('%*s', col_width, '...');
                    fprintf('%*s', col_width, X.variables{last_col});
                elseif D == 6
                    fprintf('%*s', col_width, X.variables{6});
                end
                fprintf('\n');
            
                % Display numeric values
                for i = 1:n_show
                    fprintf('    ');
                    for j = 1:D_show
                        fprintf('%*.*f', col_width, 4, X.comp(i, j));
                    end
                    if D > 6
                        fprintf('%*s', col_width, '...');
                        fprintf('%*.*f', col_width, 4, X.comp(i, last_col));
                    elseif D == 6
                        fprintf('%*.*f', col_width, 4, X.comp(i, 6));
                    end
                    fprintf('\n');
                end
            end
        else
            invalid_rows = find(errors >= tolerance);
            fprintf('✗ Some rows are NOT valid pairwise scores for D=%d\n', D);
            fprintf('  Max reconstruction error: %.2e (tolerance: %.2e)\n', ...
                error_norm, tolerance);
            fprintf('  Number of invalid rows: %d / %d\n', numel(invalid_rows), n);
            fprintf('  Row indices: %s\n', mat2str(invalid_rows'));
        end
    end
end