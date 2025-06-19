% Author: Weicen Dong <orange233@sjtu.edu.cn>
% Created Date: 2025/1/21

t1=1;a = 0.1; b = 0.1; c = 0.2; d = 0.1;
L=301;

[C,k11,w,B] = Chern_mol(t1,a,b,c,d,L);

for m=1:L+1
    for n=1:L+1
        if abs(real(shiftdim(w(2,m,n)))-real(shiftdim(w(1,m,n))))<10^(-10),rd(m,n)=1;else,rd(m,n)=NaN;end
        if abs(imag(shiftdim(w(2,m,n)))-imag(shiftdim(w(1,m,n))))<10^(-10),ri(m,n)=1;else,ri(m,n)=NaN;end
    end
end

[x,y]=meshgrid(k11,k11);

xte=-0.15; yte=1;
figure('Color','white')
subplot(2,2,1)

surf(k11/2/pi,k11/2/pi,real(shiftdim(w(1,:,:)))),hold on
surf(k11/2/pi,k11/2/pi,real(shiftdim(w(2,:,:))))
shading interp,xlabel('k_1'),ylabel('k_2'),zlabel('Re(E)')
rd=rd.*real(shiftdim(w(1,:,:)));
plot3(x(:)/2/pi,y(:)/2/pi,rd(:),'r.')%把矩阵变成向量
zlim([-4,5]),view(-20,10),hold off
 text(0.45, 1.15, '$a=b$', 'Units', 'normalized', 'FontSize', 16, 'HorizontalAlignment', 'center','Color',"#D95319",  'Interpreter', 'latex');


subplot(2,2,3)
surf(k11/2/pi,k11/2/pi,imag(shiftdim(w(1,:,:)))),hold on
surf(k11/2/pi,k11/2/pi,imag(shiftdim(w(2,:,:))))
shading interp,xlabel('k_1'),ylabel('k_2'),zlabel('Im(E)')
zlim([-1,1]),view(-20,10),hold off
zticks([-1 0 1])


t1=1;a = 0.1; b = 0.1; c = 0.2; d = 0.2;
L=301;
[C,k11,w,B] = Chern_mol(t1,a,b,c,d,L);

for m=1:L+1
    for n=1:L+1
        if abs(real(shiftdim(w(2,m,n)))-real(shiftdim(w(1,m,n))))<10^(-10),rd(m,n)=1;else,rd(m,n)=NaN;end
        if abs(imag(shiftdim(w(2,m,n)))-imag(shiftdim(w(1,m,n))))<10^(-10),ri(m,n)=1;else,ri(m,n)=NaN;end
    end
end

[x,y]=meshgrid(k11,k11);

subplot(2,2,2)
surf(k11/2/pi,k11/2/pi,real(shiftdim(w(1,:,:)))),hold on
surf(k11/2/pi,k11/2/pi,real(shiftdim(w(2,:,:))))
shading interp,xlabel('k_1'),ylabel('k_2'),zlabel('Re(E)')
rd=rd.*real(shiftdim(w(1,:,:)));
plot3(x(:)/2/pi,y(:)/2/pi,rd(:),'r.')%把矩阵变成向量
zlim([-4,5]),view(-20,10),hold off
text(0.45, 1.15, '$a\neq b$', 'Units', 'normalized', 'FontSize', 16, 'HorizontalAlignment', 'center','Color',	"#0072BD",  'Interpreter', 'latex');


subplot(2,2,4)
surf(k11/2/pi,k11/2/pi,imag(shiftdim(w(1,:,:)))),hold on
surf(k11/2/pi,k11/2/pi,imag(shiftdim(w(2,:,:))))
shading interp,xlabel('k_1'),ylabel('k_2'),zlabel('Im(E)')
zlim([-1,1]),view(-20,10),hold off
zticks([-1 0 1])


