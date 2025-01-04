function [dG,d2G] = WCA(r2,E_int,sigma)
  E = E_int;
  sigma2 = sigma^2;
  sz = size(r2);
  r2 = r2/sigma2;
%   ind = r2 < 2^(1/3);
  ind = (r2 > 0) & (r2 < 2^(1/3));
  r2 = r2(ind);
  r2i = 1./r2;
  r6i = r2i.^3;
  r12i = r6i.^2;
  dG = zeros(sz);
  dG(ind) = (24 * E * r2i .* (-r12i + r6i/2)) / sigma2;
  if nargout>1
    d2G = zeros(sz);
    d2G(ind) = (24 * E * r2i.^2 .* (7*r12i - 2*r6i)) / sigma2^2;
  end
end
