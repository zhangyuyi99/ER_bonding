%%% This is a simulation of molecule movement at the ER–mitochondria contact site 

% Simulation parameters

dt       = 3e-5;  % time step
nSteps   = 1000;  % total number of time steps
nDim     = 1;     % dimension 

L_ER     = 10;  % length of ER (um)
L_contact= 2;  % length of contact site (um)
nMol     = 100;  % number of molecules

posMol   = zeros(nMol,nDim);  % positions of molecules
isBound  = zeros(nMol,nDim);  % array/logical indicating if each molecule is bound (true/false), 
                              % all molecules are unbound by default at the start 

D        = 1;  % diffusion coefficient
kb       = 0.5;  % binding (association) probability or rate
ku       = 0.3;  % unbinding (dissociation) probability or rate

% (Any additional parameters for interactions, e.g. potential or pairwise force)
E_int = 1;  % interaction energy scale

% Options for initial distribution
initDist = 'uniform';  % or 'random', etc.

% Visualization / plotting toggles
doPlot   = true;  % or false to disable intermediate plotting
plotStep = 100;   % how often (every N steps) to update the plot


% Initialize the position 
posMol = initPositions(nMol, L_contact, initDist, nDim);


% Running simulation
for st = 1:step

end
