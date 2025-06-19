function [E,H,V] = H_zig_mol(k,t1,a,b,c,d,n,m)
% Author: Weicen Dong <orange233@sjtu.edu.cn>
% Created Date: 2025/1/21

%2n subs for a strip of honeycomb plane (zigzag OBC)
f31=2*cos(k/2);f2=1;
bb=[];cc=[];dd=[];mm=[];
g1=(a+c+1i*(b+d))*cos(k);p1=((b-d)+1i*(c-a))*sin(k);
g2=(a+c+1i*(b+d))*cos(k/2);p2=((b-d)+1i*(c-a))*sin(-k/2);
for i=1:n-1
    bb=[bb f31 f2];cc=[cc g2+p2 g2-p2];
    dd=[dd g1+p1 g1-p1]; mm=[mm m -m];
end
bb=[bb f31];dd=[dd g1+p1 g1-p1];mm=[mm m -m];
H1=t1*(diag(bb,1)+diag(bb,-1));
H2=(diag(dd)+diag(cc,2)+diag(cc,-2));
Hm=diag(mm);

H=H1+H2+Hm;

[V,D]=eig(H);
E=diag(D);

[sortedE, idx] = sort(real(E));%按照虚数部分排序
E = E(idx);
V = V(:,idx);
end