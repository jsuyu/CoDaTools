function is_SBP = check_sbp(SBP, tolerance)
% is_sbp checks if matrix is a valid Sequential Binary Partition (SBP)
%
% Verifies that SBP is a valid with correct dimensions,
% values, partitions, and resulting orthogonal basis
%
% Syntax:
%   is_SBP = check_sbp(SBP)
%   is_SBP = check_sbp(SBP, tolerance)
%
% Input:
%  SBP - table or matrix (D x (D-1))
%         Should contain only -1, 0, and +1
%   tolerance - (optional) Orthogonality tolerance (default: 1e-12)
%
% Outputs:
%   is_SBP - TRUE if SBP is a valid SBP matrix
%
% Checks performed:
%   1. Dimensions: SBP must have D rows and D-1 columns
%   2. Values: only contains -1, 0, +1
%   3. Partitions: each columns has at least one +1 and one -1
%   4. Orthogonality: all columns are mutually orthogonal (dot product ≈ 0)
%
% Example:
%   SBP = [ 1,  1,  0;
%           1, -1,  0;
%          -1,  0,  1;
%          -1,  0, -1];
%   is_SBP = check_sbp(SBP);
%
% Author: Jordi
% Date: 2025
% Updated: 2025

    % Validate tolerance and display
    if nargin < 2 || isempty(tolerance)
        tolerance = 1e-12;
    end
    
    if ~isnumeric(tolerance) || ~isscalar(tolerance) || tolerance < 0
        error('check_sbp:InvalidTolerance', 'Tolerance must be a non-negative scalar.');
    end
    
    % Validate SBP
    if istable(SBP)
        SBP = table2array(SBP);
    end

    if ~isnumeric(SBP)
        error('check_sbp:InvalidInput', 'Input must be numeric');
    end
    
    % Get dimensions
    [D, K] = size(SBP);
    
    % CHECK 1: Dimensions (K = D-1)
    if K ~= D - 1
        error('check_sbp:InvalidDimensions', ...
              'SBP must have D rows and D-1 columns. Found %d × %d.', D, K);
    end
        
    % CHECK 2: Values (only -1, 0, +1)
    unique_vals = unique(SBP(:));
    if ~all(ismember(unique_vals, [-1, 0, 1]))
        error('check_sbp:InvalidValues', ...
              'SBP must contain only -1, 0, or +1. Found: %s', mat2str(unique_vals'));
    end
    
    % CHECK 3: Each partition valid (has both + and -)
    invalid_partitions = [];
    for i = 1:K
        n_pos = sum(SBP(:, i) == 1);
        n_neg = sum(SBP(:, i) == -1);
        
        if n_pos == 0 || n_neg == 0
            invalid_partitions = [invalid_partitions, i];
        end
    end
    
    if ~isempty(invalid_partitions)
        error('check_sbp:InvalidPartitions', ...
              'Each SBP column must have at least one +1 and one -1. Invalid columns: %s', ...
              mat2str(invalid_partitions));
    end
    
    % CHECK 4: Orthogonality
    nonOrthPairs = [];
    for i = 1:K-1
        for j = i+1:K
            dp = SBP(:,i)' * SBP(:,j); % dot product
            if abs(dp) > tolerance
                nonOrthPairs = [nonOrthPairs; i, j, dp]; %#ok<AGROW>
            end
        end
    end

    if ~isempty(nonOrthPairs)
        pairsStr = sprintf('(%d,%d)=%.2e ', nonOrthPairs');
        error('check_sbp:NotOrthogonal', ...
              'SBP columns not orthogonal. Non-orthogonal pairs: %s', pairsStr);
    end    
    
    % All checks passed
    is_SBP = true;
    
end