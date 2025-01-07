function posMol = initPositions(nMol, L_ER, initDist, nDim)
% INITPOSITIONS Initialize the positions of molecules in 1D or 2D.
%
%   posMol = INITPOSITIONS(nMol, L_ER, initDist, nDim)
%
%   INPUTS:
%       nMol       - Number of molecules (integer).
%       L_ER       - Length of the ER (micrometers).
%       initDist   - 'uniform' or 'random' (controls placement distribution).
%       nDim       - 1 or 2, specifying dimensionality.
%
%   OUTPUT:
%       posMol     - Positions of molecules.
%                    If nDim=1, posMol is (nMol x 1).
%                    If nDim=2, posMol is (nMol x 2).
%
%   EXAMPLES:
%       posMol = initPositions(100, 5.0, 'uniform', 1);  % 1D, uniform
%       posMol = initPositions(100, 5.0, 'random', 2);   % 2D, random

    % Validate nDim input
    if ~(nDim == 1 || nDim == 2)
        error('initPositions:InvalidDimension', ...
              'nDim must be 1 or 2, but got %d.', nDim);
    end

    % Initialize posMol based on dimension
    if nDim == 1
        % 1D: posMol is nMol x 1
        posMol = zeros(nMol, 1);
    else
        % 2D: posMol is nMol x 2
        posMol = zeros(nMol, 2);
    end

    switch lower(initDist)
        case 'uniform'
            if nDim == 1
                % Uniform in 1D: from 0 to L_ER
                posMol = linspace(0, L_ER, nMol+1).';  % Creates nMol+1 points
                posMol = posMol(1:end-1);             % Drop the final point (L_ER)
            else
                % Uniform in 2D
                % We place molecules in a regular grid from (0,0) to (L_ER, L_ER).
                % For simplicity, we can do sqrt(nMol) points in x and y each
                % (requires nMol to be a perfect square) or make an approximate grid.
                %
                % Here is one simple approach: we fill by rows/columns.
                nPoints = ceil(sqrt(nMol)); 
                xGrid   = linspace(0, L_ER, nPoints); 
                yGrid   = linspace(0, L_ER, nPoints); 
                [X, Y]  = meshgrid(xGrid, yGrid);
                % Flatten:
                gridPoints = [X(:), Y(:)];
                % If we have more grid points than needed, just truncate:
                posMol     = gridPoints(1:nMol, :);
            end
        case 'step'
            if nDim == 1
                % We want to place all particles uniformly in a 2-um region
                % centered in the middle of [0, L_ER].
                % So the left edge is (L_ER - 2)/2, the right edge is (L_ER + 2)/2.
                
                leftEdge  = (L_ER - 2)/2;        % start of the 2 um region
                rightEdge = leftEdge + 2;        % end of the 2 um region
                
                % Create nMol+1 points from leftEdge..rightEdge
                posMol = linspace(leftEdge, rightEdge, nMol).';
            else
                error('step mode is only defined for 1D in this snippet');
            end

        case 'delta'
            if nDim == 1
                % We want to place all particles in the middle of L_ER
                
                % Create nMol+1 points from leftEdge..rightEdge
                posMol = 0.5*L_ER*ones(nMol,1);
            else
                error('delta mode is only defined for 1D in this snippet');
            end
        case 'random'
            if nDim == 1
                % Random in 1D: uniform random from 0 to L_ER
                posMol = L_ER * rand(nMol, 1);
            else
                % Random in 2D: each coordinate uniform from 0 to L_ER
                posMol = L_ER * rand(nMol, 2);
            end

        otherwise
            error('initPositions:UnknownDist', ...
                  'Unknown initialization mode: %s', initDist);
    end

end
