function [iList, jList, dist, dist2] = neighborlist_1d(x, cutoff, cellSize, L_ER)
% NEIGHBORLIST_1D  Return neighbor pairs (i,j) whose separation is within cutoff in 1D.
%
%   [iList, jList, dist] = neighborlist_1d(x, cutoff, cellSize, L_ER)
%
%   INPUTS:
%       x        - (N x 1) array of particle positions in 1D
%       cutoff   - distance cutoff. Pairs with distance <= cutoff are neighbors
%       cellSize - size of each sub-cell (e.g. cutoff/2)
%       L_ER     - total length of the 1D system (for periodic boundary, if needed)
%
%   OUTPUTS:
%       iList, jList - row vectors of neighbor indices (i < j)
%       dist         - distances between these neighbors
%
%   NOTE:
%     - If you want periodic boundaries, uncomment the code that wraps positions
%       or does a shift. If not, you can skip or adapt.
%     - Make sure x is within [0, L_ER].
%     - The code is O(N) if cellSize is chosen ~ cutoff, but strictly speaking
%       still O(N) or O(N log N) depending on distribution.

    N = numel(x);
    iList = [];
    jList = [];
    dist  = [];
    dist2  = [];

    % Optionally ensure positions are within [0, L_ER], if periodic:
    % x = mod(x, L_ER);

    % 1) Number of cells
    nCell = floor(L_ER / cellSize);
    if nCell < 1
        nCell = 1;  % fallback if cellSize is too large
    end

    % 2) Recompute an exact cellSize if we want to fill the domain precisely
    cellSize = L_ER / nCell;

    % 3) Assign each particle to a cell index
    % cell index range: 0 .. nCell-1
    cellIndex = floor(x / cellSize);
    % cell index range: 1 .. nCell, +1 for MATLAB cell indexing
    cellIndex = cellIndex + 1; 
    % clamp it just in case of borderline rounding:
    cellIndex(cellIndex >= nCell) = nCell;
    cellIndex(cellIndex < 1)      = 1;

    % 4) Build a "cell -> list of particles" mapping
    % cellContent = cell(nCell, 1);
    % for p = 1:N
    %     cInd = cellIndex(p); 
    %     cellContent{cInd} = [cellContent{cInd}, p];
    % end
    cellContent = accumarray( cellIndex(:), (1:N)', [nCell, 1], @(pts){pts}, {[]} );

    % 5) For each cell, we only need to check neighbor candidates in:
    %    that cell, the cell to the left, and the cell to the right.
    %    (If you want 2 or more cells away, adapt accordingly.)
    for c = 1:nCell
        % The set of cells to check
        neighCells = c + [-1,0,1];
        neighCells = mod(neighCells,nCell);
        neighCells(neighCells == 0) = nCell;


        % neighCells = [c];
        % if c > 1
        %     neighCells = [neighCells, c-1];
        % end
        % if c < nCell
        %     neighCells = [neighCells, c+1];
        % end

        % Particles in current cell
        partC = cellContent{c};
        % For each neighbor cell, gather neighbor particles
        allNeighParts = [];
        for nc = 1:numel(neighCells)
            c2 = neighCells(nc);
            % particles_in_c2 = cellContent{c2};
            allNeighParts = [allNeighParts, cellContent{c2}'];
        end
        % allNeighParts = allNeighParts(:);

        % Now do pairwise checks between partC and allNeighParts
        % but we want i < j to avoid duplicates
        for idx1 = 1:numel(partC)
            iP = partC(idx1);
            for idx2 = 1:numel(allNeighParts)
                jP = allNeighParts(idx2);
                if jP <= iP
                    continue;
                end

                % distance
                r = x(iP) - x(jP);

                r = r - round(r/L_ER)*L_ER;
                % so that distance is the minimal image distance

                distVal = abs(r);
                if distVal <= cutoff
                    % record neighbor
                    iList = [iList, iP];
                    jList = [jList, jP];
                    dist  = [dist, r];
                    dist2  = [dist2, r^2];
                end
            end
        end
    end
end
