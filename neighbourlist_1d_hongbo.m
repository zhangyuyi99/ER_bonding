function [i, j, r2, r] = neighborlist_1d(x, cutoff, cellsz, L)
% NEIGHBORLIST_1D  Create a neighbor list in 1D for points under periodic BC.
%
%   [i, j, r2, r] = neighborlist_1d(x, cutoff, cellsz, L)
%
%   INPUTS:
%       x        - (N x 1) vector of particle positions in [0,L)  (1D)
%       cutoff   - distance threshold for neighbor pairs
%       cellsz   - nominal cell size. The function enforces cellsz >= cutoff/2
%                  then adjusts it so that an integer # of cells fits length L
%       L        - size of the periodic domain (scalar)
%
%   OUTPUTS:
%       i, j     - row vectors listing neighbor pairs (i < j), 1-based indices
%       r2       - distance^2 between i and j
%       r        - displacement from j to i (or i to j). Here it's (x(i)-x(j)).
%
%   NOTES:
%       1) Particles i and j are considered neighbors if their 1D separation is
%          <= cutoff (with periodic wrapping).
%       2) We build a Nx-by-#cells adjacency matrix using a cell-list approach:
%           - Each particle belongs to its cell plus the two adjacent cells
%             (left & right) to handle edge crossing.
%       3) Then we quickly find candidate neighbors by sparse multiplication,
%          and filter them by an exact distance check (<= cutoff).
%
%   EXAMPLE:
%       x = 10*rand(100,1);  % random positions in [0,10)
%       cutoff = 1.0;
%       [i,j,r2,r] = neighborlist_1d(x, cutoff, cutoff/2, 10);

    % Ensure x is a column
    x = x(:);
    Nx = numel(x);

    % 1) Wrap positions into [0,L) to enforce periodic boundary
    %    (If x already in [0,L), this leaves them unchanged.)
    x = x - min(x);                % shift min to 0
    x = x - L .* floor(x ./ L);    % wrap into [0,L)

    % 2) Ensure cellsz >= cutoff/2
    cellsz = max(cellsz, cutoff/2);

    % 3) Determine how many cells fit along [0,L).
    Ncell = floor(L / cellsz);
    % If for some reason Ncell < 1, just set Ncell=1 (fallback)
    Ncell = max(Ncell, 1);

    % Recompute the actual cell size so it divides domain exactly
    cellsz = L / Ncell;

    % 4) Assign each particle to an integer cell index in [0, Ncell-1]
    ind = floor(x / cellsz);
    % Guarantee no out-of-bounds due to rounding
    ind(ind >= Ncell) = Ncell-1;
    ind(ind < 0)      = 0;

    % 5) For each particle, we consider up to 3 cells: {ind, ind-1, ind+1}, wrapped
    %    We'll store them in cind, which is Nx x 3
    %    Then build a Nx-by-(Ncell) adjacency.

    % We'll have up to 3 neighboring cells in 1D: c, c-1, c+1
    nneigh = 3;
    cind   = zeros(Nx, nneigh);   % each row => the 3 possible cell indices

    % The offsets for 1D are [-1, 0, +1]
    % We'll apply periodic wrap on cell index: (cell + Ncell) mod Ncell
    off = [-1, 0, +1];
    for p = 1:Nx
        c = ind(p);
        for n = 1:nneigh
            cshift = c + off(n);
            % wrap cell index
            cshift = mod(cshift, Ncell);
            cind(p,n) = cshift;
        end
    end

    % Convert from [0..Ncell-1] to [1..Ncell] for MATLAB indexing
    cind = cind + 1;

    % 6) Build a sparse adjacency matrix pc: Nx-by-Ncell
    %    We have Nx * 3 true entries in pc
    pind = repmat((1:Nx)', 1, nneigh);  % Nx-by-3, repeating [1..Nx] in each column
    pc   = sparse(pind(:), cind(:), true, Nx, Ncell);

    % 7) Multiply pc * pc' to see which pairs share at least one cell
    %    "neighbor" is Nx-by-Nx, set to true if two points share a cell
    neighbor = triu((pc * pc') > 0, 1);   % upper-triangular, exclude diag => i<j

    [i, j] = find(neighbor);

    % 8) Compute 1D distances with minimal-image for periodic BC
    dx = x(i) - x(j);
    dx = dx - L * round(dx / L);
    dr2 = dx.^2;
    roi = dr2 <= cutoff^2;   % index of pairs truly within cutoff

    % filter out pairs beyond cutoff
    i = i(roi);
    j = j(roi);

    % Return optional outputs
    if nargout > 2
        r2 = dr2(roi);
    end
    if nargout > 3
        r  = dx(roi);
    end
end
