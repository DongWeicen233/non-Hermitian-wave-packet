% Author: Weicen Dong <orange233@sjtu.edu.cn>
% Created Date: 2025/1/21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gca=figure('Color', 'w');

xte=-0.1; yte=1.1;
n=20;
m=1801;kk=linspace(0,2*pi,m);
u=1;t1=1;a=0.1;b=0.1;c=0.2;d=-0.1;m_num=0;
for i=1:m
    k=kk(i);
   [E,H,V] = H_zig_mol(k,t1,a,b,c,d,n,m_num);HH(:,:,i)=H;ek(:,i)=E;U(:,:,i)=V;
end
for j=1:m % different k
    for i=1:n*2%different energy
        p(i,j)=0;l=U(:,i,j);
        for o=1:2*n,p(i,j)=p(i,j)+(abs(l(o))^2)*(-1+2*(o-1)/(2*n-1));end
    end
end
subplot(2,6,2),colormap jet
l1=1;l2=2*n;
for i=l1:l2
scatter(kk/2/pi, real(ek(i,:)), 8, p(i,:), 'filled'),hold on
patch([-1.1 NaN],[1.1 NaN],[-1.1 NaN],'EdgeColor','interp')
patch([-1.1 NaN],[1.1 NaN],[1.1 NaN],'EdgeColor','interp')
end
for i=n:n+1
    scatter(kk/2/pi, real(ek(i,:)), 8, p(i,:), 'filled'),hold on
    patch([-1.1 NaN],[1.1 NaN],[-1.1 NaN],'EdgeColor','interp')
    patch([-1.1 NaN],[1.1 NaN],[1.1 NaN],'EdgeColor','interp')
end
xlim([kk(1)/2/pi kk(end)/2/pi]),ylim([-4 4])
xlabel('k'),ylabel('Re(E_k)'),hold off
text(0.5,1.14 ,'$m<n, a>b$','Units', 'normalized','HorizontalAlignment', 'center', 'FontSize', 16,'Color','r',  'Interpreter', 'latex');
text(xte, 1.14, '(b)', 'Units', 'normalized', 'FontSize', 12, 'HorizontalAlignment', 'center');

subplot(2,6,8),colormap jet
for i=l1:l2
scatter(kk/2/pi, imag(ek(i,:)), 8, p(i,:), 'filled'),hold on
patch([-1.1 NaN],[1.1 NaN],[-1.1 NaN],'EdgeColor','interp')
patch([-1.1 NaN],[1.1 NaN],[1.1 NaN],'EdgeColor','interp')
end
for i=n:n+1
    scatter(kk/2/pi, imag(ek(i,:)), 8, p(i,:), 'filled'),hold on
    patch([-1.1 NaN],[1.1 NaN],[-1.1 NaN],'EdgeColor','interp')
    patch([-1.1 NaN],[1.1 NaN],[1.1 NaN],'EdgeColor','interp')
end
xlim([kk(1)/2/pi kk(end)/2/pi])%,ylim([-0.5 0.5])
xlabel('k'),ylabel('Im(E_k)')


a=0.1;b=-0.1;c=0.2;d=0.1;m_num=0;
for i=1:m
    k=kk(i);
   [E,H,V] = H_zig_mol(k,t1,a,b,c,d,n,m_num);HH(:,:,i)=H;ek(:,i)=E;U(:,:,i)=V;
end
for j=1:m % different k
    for i=1:n*2%different energy
        p(i,j)=0;l=U(:,i,j);
        for o=1:2*n,p(i,j)=p(i,j)+(abs(l(o))^2)*(-1+2*(o-1)/(2*n-1));end
    end
end
subplot(2,6,4),colormap jet
l1=1;l2=2*n;
for i=l1:l2
scatter(kk/2/pi, real(ek(i,:)), 8, p(i,:), 'filled'),hold on
patch([-1.1 NaN],[1.1 NaN],[-1.1 NaN],'EdgeColor','interp')
patch([-1.1 NaN],[1.1 NaN],[1.1 NaN],'EdgeColor','interp')
end
for i=n:n+1
    scatter(kk/2/pi, real(ek(i,:)), 8, p(i,:), 'filled'),hold on
    patch([-1.1 NaN],[1.1 NaN],[-1.1 NaN],'EdgeColor','interp')
    patch([-1.1 NaN],[1.1 NaN],[1.1 NaN],'EdgeColor','interp')
end
xlim([kk(1)/2/pi kk(end)/2/pi]),ylim([-4 4])
xlabel('k'),ylabel('Re(E_k)'),hold off
text(0.5,1.14 ,'$m<n, a<b$','Units', 'normalized','HorizontalAlignment', 'center', 'FontSize', 16,'Color','r',  'Interpreter', 'latex');
text(xte, 1.14, '(d)', 'Units', 'normalized', 'FontSize', 12, 'HorizontalAlignment', 'center');

subplot(2,6,10),colormap jet
for i=l1:l2
scatter(kk/2/pi, imag(ek(i,:)), 8, p(i,:), 'filled'),hold on
patch([-1.1 NaN],[1.1 NaN],[-1.1 NaN],'EdgeColor','interp')
patch([-1.1 NaN],[1.1 NaN],[1.1 NaN],'EdgeColor','interp')
end
for i=n:n+1
    scatter(kk/2/pi, imag(ek(i,:)), 8, p(i,:), 'filled'),hold on
    patch([-1.1 NaN],[1.1 NaN],[-1.1 NaN],'EdgeColor','interp')
    patch([-1.1 NaN],[1.1 NaN],[1.1 NaN],'EdgeColor','interp')
end
xlim([kk(1)/2/pi kk(end)/2/pi])%,ylim([-0.5 0.5])
xlabel('k'),ylabel('Im(E_k)')


a=0.2;b=-0.1;c=0.1;d=0.1;m_num=0;
for i=1:m
    k=kk(i);
   [E,H,V] = H_zig_mol(k,t1,a,b,c,d,n,m_num);HH(:,:,i)=H;ek(:,i)=E;U(:,:,i)=V;
end
for j=1:m % different k
    for i=1:n*2%different energy
        p(i,j)=0;l=U(:,i,j);
        for o=1:2*n,p(i,j)=p(i,j)+(abs(l(o))^2)*(-1+2*(o-1)/(2*n-1));end
    end
end
subplot(2,6,5),colormap jet
l1=1;l2=2*n;
for i=l1:l2
scatter(kk/2/pi, real(ek(i,:)), 8, p(i,:), 'filled'),hold on
patch([-1.1 NaN],[1.1 NaN],[-1.1 NaN],'EdgeColor','interp')
patch([-1.1 NaN],[1.1 NaN],[1.1 NaN],'EdgeColor','interp')
end
for i=n:n+1
    scatter(kk/2/pi, real(ek(i,:)), 8, p(i,:), 'filled'),hold on
    patch([-1.1 NaN],[1.1 NaN],[-1.1 NaN],'EdgeColor','interp')
    patch([-1.1 NaN],[1.1 NaN],[1.1 NaN],'EdgeColor','interp')
end
xlim([kk(1)/2/pi kk(end)/2/pi]),ylim([-4 4])
xlabel('k'),ylabel('Re(E_k)'),hold off
text(0.5,1.14 ,'$m>n, a<b$','Units', 'normalized','HorizontalAlignment', 'center', 'FontSize', 16,'Color','r',  'Interpreter', 'latex');
text(xte, 1.14, '(e)', 'Units', 'normalized', 'FontSize', 12, 'HorizontalAlignment', 'center');
subplot(2,6,11),colormap jet
for i=l1:l2
scatter(kk/2/pi, imag(ek(i,:)), 8, p(i,:), 'filled'),hold on
patch([-1.1 NaN],[1.1 NaN],[-1.1 NaN],'EdgeColor','interp')
patch([-1.1 NaN],[1.1 NaN],[1.1 NaN],'EdgeColor','interp')
end
for i=n:n+1
    scatter(kk/2/pi, imag(ek(i,:)), 8, p(i,:), 'filled'),hold on
    patch([-1.1 NaN],[1.1 NaN],[-1.1 NaN],'EdgeColor','interp')
    patch([-1.1 NaN],[1.1 NaN],[1.1 NaN],'EdgeColor','interp')
end
xlim([kk(1)/2/pi kk(end)/2/pi])%,ylim([-0.5 0.5])
xlabel('k'),ylabel('Im(E_k)')


a=0.2;b=0.1;c=0.1;d=-0.1;m_num=0;
for i=1:m
    k=kk(i);
   [E,H,V] = H_zig_mol(k,t1,a,b,c,d,n,m_num);HH(:,:,i)=H;ek(:,i)=E;U(:,:,i)=V;
end
for j=1:m % different k
    for i=1:n*2%different energy
        p(i,j)=0;l=U(:,i,j);
        for o=1:2*n,p(i,j)=p(i,j)+(abs(l(o))^2)*(-1+2*(o-1)/(2*n-1));end
    end
end
subplot(2,6,3),colormap jet
l1=1;l2=2*n;
for i=l1:l2
scatter(kk/2/pi, real(ek(i,:)), 8, p(i,:), 'filled'),hold on
patch([-1.1 NaN],[1.1 NaN],[-1.1 NaN],'EdgeColor','interp')
patch([-1.1 NaN],[1.1 NaN],[1.1 NaN],'EdgeColor','interp')
end
for i=n:n+1
    scatter(kk/2/pi, real(ek(i,:)), 8, p(i,:), 'filled'),hold on
    patch([-1.1 NaN],[1.1 NaN],[-1.1 NaN],'EdgeColor','interp')
    patch([-1.1 NaN],[1.1 NaN],[1.1 NaN],'EdgeColor','interp')
end
xlim([kk(1)/2/pi kk(end)/2/pi]),ylim([-4 4])
xlabel('k'),ylabel('Re(E_k)'),hold off
text(0.5,1.14 ,'$m>n, a>b$','Units', 'normalized','HorizontalAlignment', 'center', 'FontSize', 16,'Color','r',  'Interpreter', 'latex');
text(xte, 1.14, '(c)', 'Units', 'normalized', 'FontSize', 12, 'HorizontalAlignment', 'center');

subplot(2,6,9),colormap jet
for i=l1:l2
scatter(kk/2/pi, imag(ek(i,:)), 8, p(i,:), 'filled'),hold on
patch([-1.1 NaN],[1.1 NaN],[-1.1 NaN],'EdgeColor','interp')
patch([-1.1 NaN],[1.1 NaN],[1.1 NaN],'EdgeColor','interp')
end
for i=n:n+1
    scatter(kk/2/pi, imag(ek(i,:)), 8, p(i,:), 'filled'),hold on
    patch([-1.1 NaN],[1.1 NaN],[-1.1 NaN],'EdgeColor','interp')
    patch([-1.1 NaN],[1.1 NaN],[1.1 NaN],'EdgeColor','interp')
end
xlim([kk(1)/2/pi kk(end)/2/pi])%,ylim([-0.5 0.5])
xlabel('k'),ylabel('Im(E_k)')



a=0.15;b=0.1;c=0.15;d=-0.1;m_num=0;
for i=1:m
    k=kk(i);
   [E,H,V] = H_zig_mol(k,t1,a,b,c,d,n,m_num);HH(:,:,i)=H;ek(:,i)=E;U(:,:,i)=V;
end
for j=1:m % different k
    for i=1:n*2%different energy
        p(i,j)=0;l=U(:,i,j);
        for o=1:2*n,p(i,j)=p(i,j)+(abs(l(o))^2)*(-1+2*(o-1)/(2*n-1));end
    end
end
subplot(2,6,6),colormap jet
l1=1;l2=2*n;
for i=l1:l2
scatter(kk/2/pi, real(ek(i,:)), 8, p(i,:), 'filled'),hold on
patch([-1.1 NaN],[1.1 NaN],[-1.1 NaN],'EdgeColor','interp')
patch([-1.1 NaN],[1.1 NaN],[1.1 NaN],'EdgeColor','interp')
end
for i=n:n+1
    scatter(kk/2/pi, real(ek(i,:)), 8, p(i,:), 'filled'),hold on
    patch([-1.1 NaN],[1.1 NaN],[-1.1 NaN],'EdgeColor','interp')
    patch([-1.1 NaN],[1.1 NaN],[1.1 NaN],'EdgeColor','interp')
end
xlim([kk(1)/2/pi kk(end)/2/pi]),ylim([-4 4])
xlabel('k'),ylabel('Re(E_k)'),hold off
text(0.5,1.14 ,'$m=n, a>b$','Units', 'normalized','HorizontalAlignment', 'center', 'FontSize', 16,'Color','r',  'Interpreter', 'latex');
text(xte, 1.14, '(f)', 'Units', 'normalized', 'FontSize', 12, 'HorizontalAlignment', 'center');

subplot(2,6,12),colormap jet
for i=l1:l2
scatter(kk/2/pi, imag(ek(i,:)), 8, p(i,:), 'filled'),hold on
patch([-1.1 NaN],[1.1 NaN],[-1.1 NaN],'EdgeColor','interp')
patch([-1.1 NaN],[1.1 NaN],[1.1 NaN],'EdgeColor','interp')
end
for i=n:n+1
    scatter(kk/2/pi, imag(ek(i,:)), 8, p(i,:), 'filled'),hold on
    patch([-1.1 NaN],[1.1 NaN],[-1.1 NaN],'EdgeColor','interp')
    patch([-1.1 NaN],[1.1 NaN],[1.1 NaN],'EdgeColor','interp')
end
xlim([kk(1)/2/pi kk(end)/2/pi])%,ylim([-0.5 0.5])
xlabel('k'),ylabel('Im(E_k)')

cb = colorbar;cb.Ticks = [-1,0, 1]; cb.TickLabels = {'upper edge','bulk', 'lower edge'}; % Set the tick labels
cb.Position=[0.92,0.2,0.01,0.6];cb.FontSize=12; cb.YDir = 'reverse';% Position and size of the colorbar

% %zigzag
subplot(1,6,1)
Nx=30;Ny=8;a_cell=1;
ratio=5;%画晶格位置最后的放缩比例，使得画面不要超出figure
   hex_x = [0, 0.5*sqrt(3), 0.5*sqrt(3), 0, -0.5*sqrt(3), -0.5*sqrt(3),0] * a_cell;
    hex_y = [-1, -0.5, 0.5, 1, 0.5, -0.5,-1] * a_cell;
    hold on;axis equal;
    for i=1:Ny-1
        for j=1:(Nx-1)/2
        plot((hex_x+(j-1+(i-1)/2)*sqrt(3))/ratio, (hex_y-(i-1)*1.5)/ratio, 'LineWidth', 0.5,'Color',[0.3 0.3 0.3]); 
        end
    end
    axis off,xlim([2.5 4.1])

x_s=(2.94449+3.11769)/2;
plot([x_s x_s],[-2.05 0.25],'b--','LineWidth',1.5)
plot([x_s+sqrt(3)/ratio x_s+sqrt(3)/ratio],[-2.05 0.25],'b--','LineWidth',1.5)
plot([x_s x_s+sqrt(3)/ratio],[-2.05 -2.05],'b--','LineWidth',1.5)
plot([x_s+sqrt(3)/ratio x_s],[0.25 0.25],'b--','LineWidth',1.5)
text(x_s+sqrt(3)/ratio/2,0.33,'unit cell', 'FontSize', 12,'HorizontalAlignment', 'center')
text(xte+0.1, yte-0.2, '(a)', 'Units', 'normalized', 'FontSize', 12, 'HorizontalAlignment', 'center');

text(x_s+sqrt(3)/ratio/2,-2.15, '$\mathbf{a}_x=(1,0)$',  'Interpreter', 'latex','HorizontalAlignment', 'center', 'FontSize', 12);

