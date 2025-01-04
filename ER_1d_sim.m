%%% This is a simulation of molecule movement at the ER–mitochondria contact site 

% Simulation parameters


param.steps    = 1000;  % total number of time steps
param.ndims    = 1;     % dimension 

param.L_ER     = 10;  % length of ER (um)
param.L_contact= 2;  % length of contact site (um)
param.sigma    = 1e-1; % size of the binding molecules 
param.N        = 1000;  % number of molecules

param.Dt       = param.sigma^2 / 3; %translational diffusivity
param.kb       = 0.5;  % binding (association) probability or rate
param.ku       = 0.3;  % unbinding (dissociation) probability or rate

param.dt       = 1e-5 * param.sigma^2 / param.Dt;  % time step
noise_strength = sqrt(2 * param.Dt * param.dt);

kBT = 1;  % Thermal energy

% (Any additional parameters for interactions, e.g. potential or pairwise force)
param.E_int = 1;  % interaction energy scale

% Options for initial distribution
param.initDist = 'uniform';  % or 'random', etc.

% Visualization / plotting toggles
param.doPlot   = true;  % or false to disable intermediate plotting
param.plotStep = 100;   % how often (every N steps) to update the plot


x   = zeros(param.N,param.ndims);  % positions of molecules
isBound  = zeros(param.N,param.ndims);  % array/logical indicating if each molecule is bound (true/false), 
                              % all molecules are unbound by default at the start 

% Initialize the position 
x = initPositions(param.N, param.L_ER, param.initDist, param.ndims);
% Preallocate array to store positions over time
xout = zeros(param.N, param.ndims, param.steps);
xout(:,:,1) = x;

% WCA cutoff: typically r_c = 2^(1/6) * sigma
rCut    = 2^(1/6) * param.sigma; 
cellSize = rCut/2;

N = param.N;
% Running simulation
for st = 2:param.steps
    % Shift x to within the box (periodic boundary conditions)
    x = x - param.L_ER  * floor(x / param.L_ER);  
    [neigh_i,neigh_j,r,r2] = neighbourlist_1d(x, rCut, cellSize, param.L_ER);
    
    % Prepare displacement
    dx = zeros(size(x));

    dG = feval("WCA",r2,param.E_int,param.sigma);

    dGdri = 2*r.*dG;
    F = accumarray(neigh_i(:),-dGdri(:),[N,1]) + accumarray(neigh_j(:),dGdri(:),[N,1]);

    dx = dx + (param.Dt/kBT) * F * param.dt; 
    % Add noise
    % dx = dx + noise_strength .* randn(param.N, param.ndims);
    % Update particle positions
    x = x + dx;
    xout(:,:,st) = x;
end

save('test_result.mat', 'xout', 'param');