% Author: Weicen Dong <orange233@sjtu.edu.cn>
% Created Date: 2025/1/21

clear all

tic
  E_fun = @(kx,ky) 0.3*kx+0.1*ky+1i*(0.1*kx+0.2*ky);
 dE1_dkx_func =@(kx,ky) 0.3+0.1i;
 dE1_dky_func =@(kx,ky) 0.1+0.2i;

% 2D grid
L=300;ky_vec=linspace(-1.5,2,L);kx_vec=linspace(-1.5,2.5,L);
ky_vals=ky_vec;kx_vals=kx_vec;

[KX,KY] = meshgrid(kx_vec,ky_vec);
z = E_fun(KX,KY);
figure,subplot(1,2,1),surf(KX,KY,real(z)); 
xlabel('kx'),ylabel('ky'),title('Re(E)'),shading interp
subplot(1,2,2),surf(KX,KY,imag(z));
xlabel('kx'),ylabel('ky'),title('Im(E)'),shading interp

%F=[0 0];
F=[0.2 0];
B=1+1i*0.5;
%B=0;

kx_bar=0;ky_bar=0;Delta_x=0.5;Delta_y=0.5;%因为ky的单位长度比kx的短sqrt(3)倍，这里我们要让k空间上两个方向delta一样。
x_k_0=@(kx,ky) exp(-((ky-ky_bar)/Delta_y).^2/2).*exp(-((kx-kx_bar)/Delta_x).^2/2)/2*pi/Delta_x/Delta_y;
Et=@(kx,ky,t) integral(@(tau) E_fun(kx+F(1)*tau, ky+F(2)*tau), 0, t);
gamma=@(kx,ky,t) -B*integral(@(tau) ky+F(2)*tau, 0, t)*F(1);
x_k_t=@(kx,ky,t) exp(-((ky-ky_bar-F(2)*t)/Delta_y).^2/2).*exp(-((kx-kx_bar-F(1)*t)/Delta_x).^2/2)/2*pi/Delta_x/Delta_y.*exp(-1i*arrayfun(@(kx,ky)Et(kx,ky,t)-gamma(kx,ky,t),kx,ky));
x_start=0;y_start=0;r_max_the(1,1)=x_start;r_max_the(1,2)=y_start;

nt=16;t_vec=linspace(0,5,nt);nt=length(t_vec);dt=t_vec(end)-t_vec(end-1);%时间再多的话k空间位置也不会动，倾向于delta function，实空间就倾向于平面波，求最大值就不准了
k_max_the=zeros(nt,2);
k_max_the(1,:)=[kx_bar ky_bar];
gaussian2D = @(p, x, y) p(1) * exp( -((x - p(2))/p(4)).^2/2  -((y - p(3))/p(5)).^2/2 );


for i=1:nt
    t=t_vec(i);
disp(i)
if i==1
    for m=1:L
        for n=1:L
    X_K_T(m,n,i)=x_k_0(KX(m,n), KY(m,n));
        end
    end
else
    for m=1:L
        for n=1:L
    X_K_T(m,n,i)=x_k_t(KX(m,n), KY(m,n), t);
        end
    end
end
     [max_val, max_index] = max(abs(X_K_T(:,:,i)), [], 'all', 'linear'); % 计算最大值
    [max_m, max_n] = ind2sub(size(X_K_T(:,:,i)), max_index); % 将线性索引转换为行列索引
    k_max(i, :) = [KX(max_m, max_n), KY(max_m,max_n)]; % 存储对应的 (kx, ky)

        % 初始猜测：峰值位置和宽度
        p0 = [max_val, k_max(i, 1),k_max(i, 2), 0.5,0.5]; % [幅值, k中心, k宽度
        % 准备数据进行拟合
        X_data = abs(X_K_T(:,:,i)); % 数据
        X_data = X_data(:); % 将二维数据展开为一维向量
        % 使用 lsqcurvefit 进行拟合
        options = optimset('Display', 'off'); % 不显示输出
        p_fit = lsqcurvefit(@(p, xy) gaussian2D(p, xy(:,1), xy(:,2)), p0, [KX(:), KY(:)], X_data, [], [], options);

        % 得到拟合结果的宽度 Delta
        Deltax(i) = p_fit(4);Deltay(i) = p_fit(5);
        k_max2(i, :) = [p_fit(2) p_fit(3)];

end

v = VideoWriter(sprintf('2D_k.mp4'), 'MPEG-4');
rate=5; v.FrameRate = rate;open(v);
figure('Color','white')
for i=1:nt
    t=t_vec(i);
    pcolor(KX, KY, abs(X_K_T(:,:,i))); % 绘制x_k_t的绝对值
    shading interp; % 插值光滑
    colorbar; 
axis equal    
xlabel('kx'); ylabel('ky'); % y轴标签
    title(['Time = ' num2str(t)]); % 标题显示当前时间
    frame = getframe(gcf);writeVideo(v,frame);
    pause(0.01);
end
close(v);




phi_r_t=@(x,y,j) trapz(ky_vals,trapz(kx_vals,X_K_T(:,:,j).*exp(1i*KX.*(x-x_start)+1i*KY.*(y-y_start)),2),1);
N_kx=40;
[x_vec,y_vec]=meshgrid(linspace(-4.5,6,N_kx),linspace(-4.5,6,N_kx));
for i=1:nt
    disp(i)
    for m=1:N_kx
        for n=1:N_kx
            phi_temp(m,n)=phi_r_t(x_vec(m,n),y_vec(m,n),i);%m:ky不同，n:kx不同
        end
    end
phi_r(:,:,i)=phi_temp;
end

v = VideoWriter(sprintf('video_phi_r_t_try.mp4'), 'MPEG-4');
rate=5; v.FrameRate = rate;open(v);
figure('Color','white')
for i=1:nt
pcolor(x_vec, y_vec, abs(phi_r(:,:,i))); % 绘制phi_r_t的绝对值
shading interp; % 插值光滑
colorbar; axis equal,xlim([min(min(x_vec)) max(max(x_vec))]),ylim([min(min(y_vec)) max(max(y_vec))])
xlabel('x'); ylabel('y'); % y轴标签
title(['Time = ' num2str(t_vec(i))]); % 标题显示当前时间
frame = getframe(gcf);writeVideo(v,frame);
pause(0.01);
[max_val, max_index] = max(abs(phi_r(:,:,i)), [], 'all', 'linear'); % 计算最大值
    [max_m, max_n] = ind2sub(size(phi_r(:,:,i)), max_index); % 将线性索引转换为行列索引
    r_max(i, :) = [x_vec(max_m, max_n), y_vec(max_m,max_n)]; 
end
close(v);

for i=2:nt
        kx=k_max_the(i-1,1);      ky=k_max_the(i-1,2);
        dt=t_vec(i)-t_vec(i-1);
        Delta_x=Deltax(i);
        Delta_y=Deltay(i); 
        kx_temp=k_max_the(i-1,1)+Delta_x^2*dt*imag(dE1_dkx_func(kx,ky))+F(1)*dt;
        ky_temp=k_max_the(i-1,2)+Delta_y^2*dt*imag(dE1_dky_func(kx,ky))+F(2)*dt;
    k_max_the(i,1)=k_max_the(i-1,1)+Delta_x^2*dt*imag(dE1_dkx_func(kx,ky)+dE1_dkx_func(kx_temp,ky_temp))/2+F(1)*dt;
    k_max_the(i,2)=k_max_the(i-1,2)+Delta_y^2*dt*imag(dE1_dky_func(kx,ky)+dE1_dky_func(kx_temp,ky_temp))/2+F(2)*dt+Delta_y^2*dt*F(1)*imag(B);
    r_max_the(i,1)=r_max_the(i-1,1)+dt*real(dE1_dkx_func(kx,ky)+dE1_dkx_func(kx_temp,ky_temp))/2;
    r_max_the(i,2)=r_max_the(i-1,2)+dt*real(dE1_dky_func(kx,ky)+dE1_dky_func(kx_temp,ky_temp))/2+dt*F(1)*real(B); 
end

figure('Color', 'w'),sgtitle('maximum of k and r')
subplot(2,2,1),plot(k_max(:,1),k_max(:,2),'o'),hold on
plot(k_max_the(:,1),k_max_the(:,2),'*')
xlabel('kx'),ylabel('ky')
subplot(2,2,2),plot(r_max(:,1),r_max(:,2),'o'),hold on
plot(r_max_the(:,1),r_max_the(:,2),'*')
legend('simulation','theory')
xlabel('x'),ylabel('y')
subplot(2,2,3),plot(t_vec,real(E_fun(k_max(:,1),k_max(:,2))))
xlabel('t'),ylabel('Re E')
subplot(2,2,4),plot(t_vec,imag(E_fun(k_max(:,1),k_max(:,2))))
xlabel('t'),ylabel('Im E')

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load('withouB.mat')
k1=k_max;
k2=k_max_the;
r1=r_max;
r2=r_max_the;
load('withB.mat')
k1B=k_max;
k2B=k_max_the;
r1B=r_max;
r2B=r_max_the;
figure('Color','white');
subplot(1,2,1),plot(k1(:,1),k1(:,2),'ob','LineWidth',1),hold on
plot(k2(:,1),k2(:,2),'Color','b','LineWidth',1)
plot(k1B(:,1),k1B(:,2),'or','LineWidth',1),hold on
plot(k2B(:,1),k2B(:,2),'Color','r','LineWidth',1)
title('k space')%,axis equal
text(-0.18, 1, '(a)', 'Units', 'normalized', 'HorizontalAlignment', 'center','FontSize',16);
xlabel('k_x'),ylabel('k_y'),set(gca,'FontSize',14)
subplot(1,2,2),plot(r1(:,1),r1(:,2),'ob','LineWidth',1),hold on
plot(r2(:,1),r2(:,2),'Color','b','LineWidth',1)
plot(r1B(:,1),r1B(:,2),'or','LineWidth',1),hold on
plot(r2B(:,1),r2B(:,2),'Color','r','LineWidth',1)
text(-0.18, 1, '(b)', 'Units', 'normalized', 'HorizontalAlignment', 'center','FontSize',16);
xlabel('x'),ylabel('y'),title('real space')
set(gca,'FontSize',14)


load('withouB.mat')
figure('Color','white')
subplot(4,2,1),plot(k_max(:,1),k_max(:,2),'ob','LineWidth',1),hold on
plot(k_max_the(:,1),k_max_the(:,2),'Color','b','LineWidth',1)
text(0.45, 1.15, 'B=0', 'Units', 'normalized', 'FontSize', 12, 'HorizontalAlignment', 'center','Color',	"#0072BD");
xlabel('kx'),ylabel('ky')
subplot(4,2,3),plot(r_max(:,1),r_max(:,2),'or','LineWidth',1),hold on
plot(r_max_the(:,1),r_max_the(:,2),'Color','r','LineWidth',1)
xlabel('x'),ylabel('y')

load('withB.mat')
subplot(4,2,2),plot(k_max(:,1),k_max(:,2),'ob','LineWidth',1),hold on
plot(k_max_the(:,1),k_max_the(:,2),'Color','b','LineWidth',1)
xlabel('kx'),ylabel('ky')
text(0.45, 1.15, 'B=1+1i*0.5', 'Units', 'normalized', 'FontSize', 12, 'HorizontalAlignment', 'center','Color',"#EDB120");
subplot(4,2,4),plot(r_max(:,1),r_max(:,2),'or','LineWidth',1),hold on
plot(r_max_the(:,1),r_max_the(:,2),'Color','r','LineWidth',1)
xlabel('x'),ylabel('y')

load('my_red.mat'),load('my_dark.mat');
load('withoutB_figure.mat')
nx=0.5;ny=0.9;k_xmmin=-1;k_xmmax=1.7;k_ymmin=-1;k_ymmax=1.7;
xmin=min(min(x_vec))-0.5;xmax=max(max(x_vec))+0.5;ymin=min(min(y_vec))-0.5;ymax=max(max(y_vec))+0.5;
h1=subplot(4,3,7);
i=1;
pcolor(KX, KY, abs(X_K_T(:,:,i))); colormap(h1,my_dark),axis equal ,  shading interp;
text(nx,ny, 't=0', 'Units', 'normalized', 'FontSize', 9, 'HorizontalAlignment', 'center');
xlim([k_xmmin k_xmmax]),ylim([k_ymmin k_ymmax]);
xlabel('kx'); ylabel('ky'); % y轴标签
h2=subplot(4,3,8);
i=nt;
pcolor(KX, KY, abs(X_K_T(:,:,i))); colormap(h2,my_dark),axis equal,  shading interp;
text(nx,ny, 't=5 (B=0)', 'Units', 'normalized', 'FontSize', 9, 'HorizontalAlignment', 'center');
xlim([k_xmmin k_xmmax]),ylim([k_ymmin k_ymmax]);
xlabel('kx'); ylabel('ky'); % y轴标签
h1=subplot(4,3,10);
i=1;
pcolor(x_vec, y_vec, abs(phi_r(:,:,i))); colormap(h1,my_red),axis equal ,  shading interp;
text(nx,ny, 't=0', 'Units', 'normalized', 'FontSize', 9, 'HorizontalAlignment', 'center');
xlim([xmin xmax]),ylim([ymin ymax]);xlabel('x'); ylabel('y'); % y轴标签
h2=subplot(4,3,11);
i=nt;
pcolor(x_vec, y_vec, abs(phi_r(:,:,i))); colormap(h2,my_red),axis equal,  shading interp;
text(nx,ny, 't=5 (B=0)', 'Units', 'normalized', 'FontSize', 9, 'HorizontalAlignment', 'center');
xlim([xmin xmax]),ylim([ymin ymax]);xlabel('x'); ylabel('y'); % y轴标签

load('withB_figure.mat')
nx=0.5;ny=0.9;k_xmmin=-1;k_xmmax=1.7;k_ymmin=-1;k_ymmax=1.7;
xmin=-5;xmax=6.5;ymin=-5;ymax=6.5;
h2=subplot(4,3,9);
i=nt;
pcolor(KX, KY, abs(X_K_T(:,:,i))); colormap(h2,my_dark),axis equal,  shading interp;
text(nx,ny, 't=5 (B=1+1i*0.5)', 'Units', 'normalized', 'FontSize', 9, 'HorizontalAlignment', 'center');
xlim([k_xmmin k_xmmax]),ylim([k_ymmin k_ymmax]);
xlabel('kx'); ylabel('ky'); % y轴标签
h2=subplot(4,3,12);
i=nt;
pcolor(x_vec, y_vec, abs(phi_r(:,:,i))); colormap(h2,my_red),axis equal,  shading interp;
text(nx,ny, 't=5 (B=1+1i*0.5)', 'Units', 'normalized', 'FontSize', 9, 'HorizontalAlignment', 'center');
xlim([xmin xmax]),ylim([ymin ymax]);xlabel('x'); ylabel('y'); % y轴标签