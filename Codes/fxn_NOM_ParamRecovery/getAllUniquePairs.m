function pairs = getAllUniquePairs(vec)
    % getUniquePairs generates all possible unique pairs from the vector vec.
    %
    % Input:
    %   vec  - A vector of elements (numeric or cell array of strings)
    % Output:
    %   pairs - An Nx2 matrix, where N is the number of unique pairs and
    %           each row is a pair of elements from vec, without repetition.

    n = length(vec);  % Number of elements in the input vector
    numPairs = n * (n - 1) / 2;  % Total number of unique pairs
    pairs = zeros(numPairs, 2);  % Preallocate output matrix (for numeric data)

    % Counter to track the current pair
    pairIdx = 1;

    % Generate all unique pairs
    for i = 1:n-1
        for j = i+1:n
            pairs(pairIdx, :) = [vec(i), vec(j)];
            pairIdx = pairIdx + 1;
        end
    end
end
