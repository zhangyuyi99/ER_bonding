%%% This is a simulation of molecule movement at the ER–mitochondria contact site 

% Simulation parameters
etaperstep = 10;

param.steps    = 2000;  % total number of time steps
param.ndims    = 1;     % dimension 

param.L_ER     = 10;  % length of ER (um)
param.L_contact= 2;  % length of contact site (um)
param.sigma    = 5e-2; % size of the binding molecules 
% param.sigma    = 1e-1; % size of the binding molecules 
param.N        = 5000;  % number of molecules

param.Dt       = param.sigma^2 / 3; %translational diffusivity
param.kb       = 0.5;  % binding (association) probability or rate
param.ku       = 0.3;  % unbinding (dissociation) probability or rate

param.dt       = 5e-6 * param.sigma^2 / param.Dt;  % time step
noise_strength = sqrt(2 * param.Dt * param.dt);

kBT = 1;  % Thermal energy

% (Any additional parameters for interactions, e.g. potential or pairwise force)
param.pairwise_interaction = "WCA";
param.E_int = 1e-10;  % interaction energy scale
param.max_F_dx = 0.5;

% Options for initial distribution
param.initDist = 'random';  % or 'random', 'uniform', 'step', 'delta' etc.

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
isBoundout = zeros(param.N, param.ndims, param.steps);
xout(:,:,1) = x;
isBoundout(:,:,1) = isBound;

% % WCA cutoff: typically r_c = 2^(1/6) * sigma
% rCut    = 2^(1/6) * param.sigma; 
% LJ cutoff: typically r_c = 2.5 * sigma
rCut    = 2.5 * param.sigma; 

cellSize = rCut;

N = param.N;

tk = timekeeper;
tk.istart = 1;
tt = tic;

% Running simulation
for st = 2:param.steps
    % Shift x to within the box (periodic boundary conditions)
    x = x - param.L_ER  * floor(x / param.L_ER); 

    % Prepare displacement
    dx = zeros(size(x));

    % Find the neighbouring pairs 
    [neigh_i,neigh_j,r,r2] = neighbourlist_1d(x, rCut, cellSize, param.L_ER);

    dG = feval(param.pairwise_interaction,r2,param.E_int,param.sigma);

    dGdri = 2*r.*dG;
    F = accumarray(neigh_i(:),-dGdri(:),[N,1]) + accumarray(neigh_j(:),dGdri(:),[N,1]);

    dx = dx + (param.Dt/kBT) * F * param.dt; 

    dx(dx>param.max_F_dx) = param.max_F_dx;
    dx(dx<-param.max_F_dx) = -param.max_F_dx;

    % Add noise
    dx = dx + noise_strength .* randn(param.N, param.ndims);
    % Update particle positions
    x = x + dx;
    xout(:,:,st) = x;

    %% Binding and unbinding 
    % once the particle leaves the contact site, it becomes unbinded 

    % if the unbind particles are at the middle param.L_contact, there is
    % param.kb chance for them to become binded 
    
    % for all binded particles from the last time step, there is param.ku
    % chance for them to become unbind

    % store the current binding in matrix isBoundout

    % 1) Save the old state of isBound so we can decide who was bound/unbound
    isBound_old = isBound;
    
    % 2) Identify which particles are in the contact region 
    %    (assuming 1D domain [0..L_ER], contact site is [0..L_contact]).
    inContact = (x < 0.5*(param.L_ER + param.L_contact)) ...
         & (x > 0.5*(param.L_ER - param.L_contact));
    % 3) Any particle not in the contact site becomes unbound
    %    "once the particle leaves the contact site, it becomes unbinded"
    idx_out = ~inContact;       % those outside [0,L_contact)
    isBound(idx_out) = false;
    
    % 4) If the particle is unbound and now in the contact region, 
    %    it has param.kb chance to become bound
    %    We only let them bind if they were unbound in the previous step (or we can check current).
    idx_in_unbound = inContact & ~isBound_old;  
    rbind = rand(size(x));      % random [0..1] for each particle
    bindEvents = (rbind < param.kb);
    isBound(idx_in_unbound) = bindEvents(idx_in_unbound);
    
    % 5) For all particles that were bound last step, 
    %    they have param.ku chance to become unbound
    idx_bound_prev = (isBound_old == true);  
    runbind = rand(size(x)); 
    unbindEvents = (runbind < param.ku);  % e.g. 30% chance if param.ku=0.3
    % Only apply to those that *were* bound
    idx_unbind = idx_bound_prev & unbindEvents;
    isBound(idx_unbind) = false;
    
    % 6) Finally, store the current binding state in isBoundout for this time step
    isBoundout(:,:,st) = isBound;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if mod(st, etaperstep) == 0
        tk = timekeeper(tk, st, param.steps);
    end
end

toc(tt);

save('result\test_result.mat', 'xout', 'isBoundout', 'param');