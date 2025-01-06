function [dG,d2G] = LJ(r2,E_int,sigma)
  E = E_int;
  sigma2 = sigma^2;
  r2 = r2/sigma2;
  r2i = 1./r2;
  r6i = r2i.^3;
  r12i = r6i.^2;
%   dG(ind) = (24 * E * r2i .* (-r12i + r6i/2)) / sigma2;
  dG = (24 * E * r2i .* (-r12i + r6i/2)) / sigma2;
  if nargout>1
    d2G = (24 * E * r2i.^2 .* (7*r12i - 2*r6i)) / sigma2^2;
  end
end

