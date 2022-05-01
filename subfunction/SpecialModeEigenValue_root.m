function n_mode=SpecialModeEigenValue_root(r,lam,N,V,mode1)
for kk=1:length(mode1)

data = imag(V{1,mode1(kk)}(1,:))';
f = chebfun(data);
n_mode(kk)=length(roots(f));
end

