% CHEB  compute D = differentiation matrix, x = Chebyshev grid

  function [D,r] = cheb(N,a,b)  %space from [-1,1]to[a,b]
  if nargin==1 a=-1;b=1;end
  if N==0, D=0; x=1; return, end
  x = cos(pi*(0:N)/N)'; 
  r=(b-a)/2*(x+(a+b)/(b-a));
  c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
  X = repmat(x,1,N+1);
  dX = X-X';                  
  D  = (c*(1./c)')./(dX+(eye(N+1)));      % off-diagonal entries
  D  = D - diag(sum(D')); % diagonal entries
  D=2/(b-a).*D;
