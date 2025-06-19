function [E,H,V] = H2by2_mol(k1,k2,t1,a,b,c,d)
% Author: Weicen Dong <orange233@sjtu.edu.cn>
% Created Date: 2025/1/21

%NH, inversion asymmetric
f=t1*exp(1i*(k1-k2)/3)*(1+exp(1i*k2)+exp(-1i*k1));

p=(1i*(c-a)+(b-d))*(sin(k1)+sin(k2)-sin(k1+k2));
g=(a+1i*b+c+1i*d)*(cos(k1)+cos(k2)+cos(k1+k2));

H=[g+p f;conj(f) g-p];
[V,D]=eig(H);
E=diag(D);

% Sort E based on the real part of the eigenvalues
[sortedE, idx] = sort(real(E));
E = E(idx);
V = V(:,idx);


end