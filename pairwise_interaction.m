function [F,r2,info] = pairwise_interaction(x,param,particle_type,info)
  %a standalone function that gives the force due to pairwise interaction among particles
  %provide the function handle or name to the pairwise interaction function at param.substrate.pairwise.func
  %function should return d(energy)/d(r2), where r2 is the distance squared
  %provide parameters passed into function at param.substrate.pairwise.param
  %provide cutoff distance for searching neighbors at param.substrate.pairwise.cutoff (separate from param.substrate.dFdr.cutoff). Or set automatically below
  %note that cellsz for neighborlist will be automatically set at cutoff/2
  %available pairwise interaction: 
  %'WCA', for which param.substrate.pairwise.param = [E, sigma], where E is the energy coefficient. sigma is the characteristic distance. In this case cutoff will be set automatically at sigma*2^(1/3)
  %'Gaussian_truncate', see condensate_substrate
  %'soft_repulsive', one sided Hookean. param = [k, r0], where k is the Hookean strength, r0 is eqm distance, beyond the eqm distance, force is zero, within the eqm distance, it's a Hookean spring
  %'LJ' full Lennard Jones potential. param = [E, sigma], where E is the energy coefficient. sigma is the characteristic distance. User needs to set the cutoff
  d = length(param.L);
  Nx = size(x,1);
  F = zeros(Nx,d);

  sigma = param.sigma;
  E_mat = param.interaction_strength_matrix;
  switch param.substrate.pairwise.func
  case 'WCA'
    cutoff = sigma*2^(1/6);
  case 'soft_repulsive'
    cutoff = sigma;
  otherwise
    cutoff = param.cutoff;
  end
  cellsz = cutoff/2*ones(1,d);
  [neigh_i,neigh_j,r2,r] = neighborlist(x,cutoff,cellsz,param.L);
  if nargout>1
    [dG,d2G] = feval(param.substrate.pairwise.func,neigh_i,neigh_j,r2,particle_type,E_mat,sigma);
  else
    dG = feval(param.substrate.pairwise.func,neigh_i,neigh_j,r2,particle_type,E_mat,sigma);
  end
  dGdri = 2*r.*dG;
  for i = 1:d
    F(:,i) = accumarray(neigh_i,-dGdri(:,i),[Nx,1]) + accumarray(neigh_j,dGdri(:,i),[Nx,1]);
  end
%   F = F(:);
  if nargout>1
    info.r = r;
    info.dG = dG;
    info.d2G = d2G;
    info.neigh_i = neigh_i;
    info.neigh_j = neigh_j;
  end
  end
