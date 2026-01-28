function iter_idx = fxn_getIterIdx(nIter, nJob, iJob, mode)

    % This function computes the indices of iterations assigned to a specific job
    % in a parallel processing setup.
    % Inputs:
    %   nIter - Total number of iterations
    %   nJob  - Total number of jobs
    %   iJob  - Index of the current job (1-based)
    %   mode  - (Optional) Mode of assignment: 'stride' or 'block' (default: 'stride')
    % Outputs:
    %   iter_idx - Indices of iterations assigned to the current job
    
if nargin < 4, mode = 'stride'; end
switch lower(mode)
    case 'stride'
        iter_idx = iJob:nJob:nIter;
    case 'block'
        base = floor(nIter/nJob);
        rem  = mod(nIter, nJob);
        start = (iJob-1)*base + min(iJob-1, rem) + 1;
        stop  = start + base - 1 + (iJob <= rem);
        iter_idx = start:stop;
    otherwise
        error('Unknown mode: %s', mode);
end
end
