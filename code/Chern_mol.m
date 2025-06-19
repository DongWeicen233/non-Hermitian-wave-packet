function [C,k11,w,B,g,Cg,gx,Cgx] = Chern_mol(t1,a,b,c,d,L)
% Author: Weicen Dong <orange233@sjtu.edu.cn>
% Created Date: 2025/1/21

    %calculate chern number and fidenity number of a 2D topological insulator
    %also save berry curvature and quantum metric
%non Hermitian, inversion symmetric
k11=linspace(0,2*pi,L+1);k22=k11;
for m=1:1:L+1
    for n=1:1:L+1
        k1=k11(m);k2=k22(n);
        [E,H,V] = H2by2_mol(k1,k2,t1,a,b,c,d);
        valH(:,:,m,n)=H;v(:,:,m,n)=V;w(:,m,n)=E;
    end
end

dk=2*pi/L;B=zeros(L,L);g=B;gx=B;%gy=B;
for m=1:1:L
    for n=1:1:L
        R=v(:,:,m,n);Le=inv(R);
       valD1H(:,:,m,n)=(valH(:,:,m+1,n)-valH(:,:,m,n))/dk;
       valD2H(:,:,m,n)=(valH(:,:,m,n+1)-valH(:,:,m,n))/dk;
       A1f=Le*valD1H(:,:,m,n)*R;
       A2f=Le*valD2H(:,:,m,n)*R;
       B(m,n)=B(m,n)+1i*(1/((w(1,m,n)-w(2,m,n)).^2))*(A1f(1,2)*A2f(2,1)-A2f(1,2)*A1f(2,1));
       g(m,n)=g(m,n)+(1/((w(1,m,n)-w(2,m,n)).^2))*(A1f(1,2)*A2f(2,1)+A2f(1,2)*A1f(2,1))/2;
       gx(m,n)=gx(m,n)+(1/((w(1,m,n)-w(2,m,n)).^2))*(A1f(1,2)*A1f(2,1));
      % gy(m,n)=gy(m,n)+(1/((w(1,m,n)-w(2,m,n)).^2))*(A2f(1,2)*A2f(2,1));
   %w(1,:,:)的Berry curvature
    end
end
%chern number
C=-sum(sum(B(:,:)))*dk*dk*sin(pi/3)*sqrt(3)/2*4/3/2/pi*sign(w(2,1,1)-w(1,1,1));
Cg=sum(sum(g(:,:)))*dk*dk*sin(pi/3)*sqrt(3)/2*4/3/2/pi*sign(w(2,1,1)-w(1,1,1))/2/pi;%fidenity number 的计算比陈数多除了2pi
Cgx=sum(sum(gx(:,:)))*dk*dk*sin(pi/3)*sqrt(3)/2*4/3/2/pi*sign(w(2,1,1)-w(1,1,1))/2/pi;
%Cgy=-sum(sum(gy(:,:)))*dk*dk*sin(pi/3)*sqrt(3)/2*4/3/2/pi*sign(w(2,1,1)-w(1,1,1));

end