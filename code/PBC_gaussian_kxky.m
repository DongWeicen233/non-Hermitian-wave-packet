% Author: Weicen Dong <orange233@sjtu.edu.cn>
% Created Date: 2025/1/21

% PBC gaussian
% without boundary, we have wave packet on 2D plane

clear all

%系统的参数
t1 = 1;
a = 0.2; b = 0.1; c = 0.1; d = -0.1;
syms kx ky real%注意本代码中的k都是以2*pi为周期。
% 定义符号形式的 f(k1, k2), p(k1, k2), g(k1, k2)
f_sym = t1 * exp(1i * (kx - (-kx/2+sqrt(3)*ky/2)) / 3) * (1 + exp(1i * (-kx/2+sqrt(3)*ky/2)) + exp(-1i * kx));
p_sym = (1i * (c - a) + (b - d)) * (sin(kx) + sin(-kx/2+sqrt(3)*ky/2) - sin(kx + (-kx/2+sqrt(3)*ky/2)));
g_sym = (a + 1i * b + c + 1i * d) * (cos(kx) + cos(-kx/2+sqrt(3)*ky/2) + cos(kx + (-kx/2+sqrt(3)*ky/2)));
% 定义 E1(k1, k2)
E1_sym = g_sym - sqrt(p_sym^2 + f_sym * conj(f_sym)); 
E_scalar_fun = matlabFunction(E1_sym, 'Vars', [kx, ky]);
E_fun = @(kx_vals, ky_vals) arrayfun(E_scalar_fun, kx_vals, ky_vals);
% 对 E1(k1, k2) 分别对 k1 和 k2 求导
dE1_dkx = diff(E1_sym, kx);dE1_dky = diff(E1_sym, ky);
% 将符号表达式转换为数值函数
dE1_dkx_func = matlabFunction(dE1_dkx, 'Vars', [kx, ky]);
dE1_dky_func = matlabFunction(dE1_dky, 'Vars', [kx, ky]);
% 定义 H(k1, k2)
H_sym = [g_sym + p_sym, f_sym; conj(f_sym) g_sym - p_sym]; H_fun = matlabFunction(H_sym, 'Vars', [kx, ky]);
% 定义 V1(k1, k2),对应于E1的本征态 。第一个和第二个元素对应AB子晶格
V1_sym = [(p_sym - sqrt(p_sym^2 + f_sym * conj(f_sym)))/conj(f_sym); 1]; 
V1_fun_A = matlabFunction(V1_sym(1)/norm(V1_sym), 'Vars', [kx, ky]);
V1_fun_B = matlabFunction(V1_sym(2)/norm(V1_sym), 'Vars', [kx, ky]);
L1_sym = [-conj(f_sym) p_sym+sqrt(p_sym^2+f_sym*conj(f_sym))]/2/sqrt(p_sym^2+f_sym*conj(f_sym));
dV1_dkx= diff(V1_sym, kx);dV1_dky= diff(V1_sym, ky);
Ax_scalar_fun = matlabFunction(1i*L1_sym*dV1_dkx, 'Vars', [kx, ky]);
Ax_fun = @(kx_vals, ky_vals) arrayfun(Ax_scalar_fun, kx_vals, ky_vals);
Ay_scalar_fun = matlabFunction(1i*L1_sym*dV1_dky, 'Vars', [kx, ky]);
Ay_fun = @(kx_vals, ky_vals) arrayfun(Ay_scalar_fun, kx_vals, ky_vals);


% k1 k2 和 kx ky 转换
kxy_to_k1 = @(kx,ky) kx; kxy_to_k2 = @(kx,ky) -kx/2+sqrt(3)*ky/2; 
k12_to_kx = @(k1,k2) k1; k12_to_ky = @(k1,k2) k1/sqrt(3)+2*k2/sqrt(3);

L=301;
k11=linspace(0,2*pi,L);k22=k11;
for m=1:1:L
    for n=1:1:L
        k1=k11(m);k2=k22(n);
        w(m,n)=E_fun(k1,(k1+2*k2)/sqrt(3));
    end
end
[K1,K2]=meshgrid(k11,k22);
KX=K1;KY=(K1+2*K2)/sqrt(3);
figure,subplot(1,2,1),surf(KX/2/pi,KY/2/pi,real(w))
xlabel('kx'),ylabel('ky'),title('Re(E)'),shading interp,axis equal
subplot(1,2,2),surf(KX/2/pi,KY/2/pi,imag(w))
xlabel('kx'),ylabel('ky'),title('Im(E)'),shading interp,axis equal

x_min=min(min(KX));x_max=max(max(KX))/2;y_min=min(min(KY));y_max=max(max(KY))*0.7/1.7;
N_kx=200;kx_vals = linspace(x_min, x_max, N_kx); % Discretize kx with N_kx points
N_ky=N_kx;ky_vals = linspace(y_min, y_max, N_ky); % Discretize ky with N_ky points
[kx_grid, ky_grid] = meshgrid(kx_vals, ky_vals);


F=[0.05 0];
%F=[0 0];
Ny=41;Nx=81;
x_start=((Nx-1)/2)*3/5;%,((Nx-1)/2)是第一行有多少ax，x_start=0是A sub,x_start不为0是向右移动多少Asub
y_start=-round((Ny-1)/2*5/6)-1/sqrt(3);%以ax ay为单位.(Ny-1)/2是一共有多少ay，-1/3是从最上面的原子开始算，y_start是Asubs
%k 空间波包从k1=k2=1/4运行到1/3
kx_bar=2*pi/4-0.3;ky_bar=2*pi*sqrt(3)/4;Delta_x=0.6;Delta_y=Delta_x;%因为ky的单位长度比kx的短sqrt(3)倍，这里我们要让k空间上两个方向delta一样。

%initial
x_k_0=@(kx,ky,t) exp(-((ky-ky_bar-F(2)*t)/Delta_y).^2/2).*exp(-((kx-kx_bar-F(1)*t)/Delta_x).^2/2)/2*pi/Delta_x/Delta_y;
kx_lin=@(kx,t)linspace(kx,kx+F(1)*t,1+ceil(t*3));
gamma=@(kx,ky,t)0;%此处只假设了x方向的力
Et=@(kx,ky,t) trapz(linspace(0,t,ceil(t*3)),E_fun(kx+F(1)*linspace(0,t,ceil(t*3)),ky+F(2)*linspace(0,t,ceil(t*3))));
x_k_t=@(kx,ky,t) exp(-((ky-ky_bar-F(2)*t)/Delta_y).^2/2).*exp(-((kx-kx_bar-F(1)*t)/Delta_x).^2/2)/2*pi/Delta_x/Delta_y.*exp(-1i*arrayfun(@(kx,ky)Et(kx,ky,t)-gamma(kx,ky,t),kx,ky));

%！！！！！！！！！！能量factor是积分不是直接乘t！！!
nt=17;t_vec=linspace(0,8,nt);nt=length(t_vec);dt=t_vec(end)-t_vec(end-1);%时间再多的话k空间位置也不会动，倾向于delta function，实空间就倾向于平面波，求最大值就不准了
k_max_the=zeros(nt,2);r_max_the=zeros(nt,2);k_max_the(1,:)=[kx_bar ky_bar];


%高斯波包在k空间的演化



%计算
gaussian2D = @(p, x, y) p(1) * exp( -((x - p(2))/p(4)).^2/2  -((y - p(3))/p(5)).^2/2 );

for i=1:nt
    disp(i)
   t=t_vec(i);
   if i==1
    for m=1:N_ky
        for n=1:N_kx
    X_K_T(m,n,i)=x_k_0(kx_grid(m,n), ky_grid(m,n), t);%m:ky不同，n:kx不同
        end
    end
   else
    for m=1:N_ky
        for n=1:N_kx
    X_K_T(m,n,i)=x_k_t(kx_grid(m,n), ky_grid(m,n), t);%m:ky不同，n:kx不同
        end
    end
   end
    [max_val, max_index] = max(abs(X_K_T(:,:,i)), [], 'all', 'linear'); % 计算最大值
    [max_m, max_n] = ind2sub(size(X_K_T(:,:,i)), max_index); % 将线性索引转换为行列索引
    k_max(i, :) = [kx_grid(max_m, max_n), ky_grid(max_m,max_n)]; % 存储对应的 (kx, ky) 

    % 初始猜测：峰值位置和宽度
    p0 = [max_val, k_max(i, 1),k_max(i,2), 0.4,0.4]; % [幅值, kx中心, kx宽度, ky中心, ky宽度
    % 准备数据进行拟合
    X_data = abs(X_K_T(:,:,i)); % 数据
    X_data = X_data(:); % 将二维数据展开为一维向量
    % 使用 lsqcurvefit 进行拟合
    options = optimset('Display', 'off'); % 不显示输出
    p_fit = lsqcurvefit(@(p, xy) gaussian2D(p, xy(:,1), xy(:,2)), p0, [kx_grid(:), ky_grid(:)], X_data, [], [], options);
    % 得到拟合结果的宽度 Delta
      Deltax(i) = p_fit(4);Deltay(i) = p_fit(5);
    k_max2(i, :) = [p_fit(2) p_fit(3)];
end

%画图
v = VideoWriter(sprintf('video_x_k_t_try.mp4'), 'MPEG-4');
rate=5; v.FrameRate = rate;open(v);
figure('Color','white'),load('my_sky_2.mat')
for i=1:nt
        t=t_vec(i);
    pcolor(kx_grid/2/pi, ky_grid/2/pi, abs(X_K_T(:,:,i))); % 绘制x_k_t的绝对值
    shading interp; % 插值光滑
    colorbar; axis equal,xlim([min(min(KX))/2/pi max(max(KX))/2/pi]),ylim([min(min(KY))/2/pi max(max(KY))/2/pi])
    xlabel('kx'); ylabel('ky'); % y轴标签
    text(0.2,1.6,['t = ' num2str(t)],'FontSize',12)
     hold on,plot(KX(1,:)/2/pi,KY(1,:)/2/pi,'k--'),plot(KX(end,:)/2/pi,KY(end,:)/2/pi,'k--'),plot(KX(:,1)/2/pi,KY(:,1)/2/pi,'k--'),plot(KX(:,end)/2/pi,KY(:,end)/2/pi,'k--')
    hold off,colormap(my_sky_2)
    frame = getframe(gcf);writeVideo(v,frame);
    pause(0.01);
end
close(v);


%高斯波包在实空间的演化
phi_r_t=@(x,y,j) trapz(ky_vals,trapz(kx_vals,X_K_T(:,:,j).*exp(1i*kx_grid.*(x-x_start)+1i*ky_grid.*(y-y_start)),2),1);
N_kx=N_kx*3;
[x_vec,y_vec]=meshgrid(linspace(20,33,N_kx),linspace(-20,-8,N_kx));
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

[C,k11,w,B,g,Cg,gx,Cgx] = Chern_mol(t1,a,b,c,d,L);
[X,Y]=meshgrid(k11(1:end-1),k11(1:end-1));%XY是k1 k2
%B_fun = @(kx, ky) interp2(X, Y, B, kx, -kx/2+sqrt(3)*ky/2, 'spline');
B_fun = @(kx, ky) 0;

r_max_the(1,1)=x_start;r_max_the(1,2)=y_start;


%使用general的EOM
for i=1:nt-1
X_K_T_Interp =@(kx, ky) interp2(kx_grid, ky_grid, abs(X_K_T(:, :, i)),kx,ky,'spline');
% 计算 abs_X_K_T 在网格上的导数
[dX_dkx, dX_dky] = gradient(abs(X_K_T(:, :, i)), kx_grid(1,2) - kx_grid(1,1), ky_grid(2,1) - ky_grid(1,1));
dX_dkx_Interp = @(kx, ky) interp2(kx_grid, ky_grid, dX_dkx,kx,ky,'spline');
dX_dky_Interp = @(kx, ky) interp2(kx_grid, ky_grid, dX_dky,kx,ky,'spline');
  dt=t_vec(i+1)-t_vec(i);
  vx=dE1_dkx_func(k_max_the(i,1),k_max_the(i,2));
  vy=dE1_dky_func(k_max_the(i,1),k_max_the(i,2));
kx_start = 2*pi/3; % 起始点kx
ky_start = 2*pi/sqrt(3); % 起始点ky
% 定义联合目标函数，将两个目标函数的误差平方和作为目标值
 objective_combined = @(k) ...
    (dX_dkx_Interp(k(1), k(2)) - imag(vx)*dt)^2 + ...
    (dX_dky_Interp(k(1), k(2)) - imag(vy)*dt)^2;
% 使用 fminsearch 寻找满足条件的 (kx, ky)
options = optimset('Display', 'off', 'TolX', 1e-6, 'TolFun', 1e-6);
[k_sol, fval] = fminsearch(objective_combined, [k_max(i, 1), k_max(i, 2)], options);
k_max_the(i+1,1)=k_sol(1);k_max_the(i+1,2)=k_sol(2);
  vx2=dE1_dkx_func(k_max_the(i+1,1),k_max_the(i+1,2));
  vy2=dE1_dky_func(k_max_the(i+1,1),k_max_the(i+1,2));
r_max_the(i+1,1)=r_max_the(i,1)+dt*real(vx2);
r_max_the(i+1,2)=r_max_the(i,2)+dt*real(vy2);
end


figure('Color','white')
subplot(1,2,1), plot(k_max(:,1)/2/pi,k_max(:,2)/2/pi,'o','LineWidth',2),hold on
plot(k_max_the(:,1)/2/pi,k_max_the(:,2)/2/pi,'LineWidth',2)
title('k space'),xlabel('kx'),ylabel('ky')%,ylim([0.42 0.49]),xlim([0.18 0.26])
%axis off
subplot(1,2,2),plot(r_max(:,1),r_max(:,2),'o','LineWidth',2),hold on
plot(r_max_the(:,1),r_max_the(:,2),'LineWidth',2)
title('r space'),xlabel('x'),ylabel('y'),xlim([24 26.2])
%axis off
    legend('simulation','theory')
sgtitle('Maximum of Gaussian wave packet')

  VA=V1_fun_A(kx_grid,ky_grid);VB=V1_fun_B(kx_grid,ky_grid);
  phi_r_t_A=@(x,y,j) trapz(ky_vals,trapz(kx_vals,X_K_T(:,:,j).*exp(1i*kx_grid.*(x-x_start)+1i*ky_grid.*(y-y_start)).*VA,2),1);
  phi_r_t_B=@(x,y,j) trapz(ky_vals,trapz(kx_vals,X_K_T(:,:,j).*exp(1i*kx_grid.*(x-x_start)+1i*ky_grid.*(y-y_start)).*VB,2),1);



%在晶格上计算
N_tot=Nx*2+(Ny-2)*(Nx+1);a_cell=1;[save_num,marker,y_unique] = coordinate(Nx,Ny);
phi=zeros(N_tot,nt);
parfor j=1:nt
    disp(j);
    t=t_vec(j);
        temp_phi = zeros(N_tot, 1);  % Temporary variable to hold the result of the current iteration
for i=1:N_tot
    x_fi=save_num(i,2);y_fi=save_num(i,3);%x y_fi 是画图用的坐标，要转换到ax ay对应的坐标上
    x=x_fi/sqrt(3);y=y_fi/sqrt(3);
if mod(marker(i),2)==1
        temp_phi(i) = phi_r_t_A(x,y,j);
else,temp_phi(i) = phi_r_t_B(x,y,j);

end
end
    phi(:, j) = temp_phi;  % Store the result of this iteration in the corresponding column of phi
    [max_value, max_index(j)] = max(abs(temp_phi));
end

%画图
v = VideoWriter(sprintf('video_phi_r_t_lattice.mp4'), 'MPEG-4');
rate=5; v.FrameRate = rate;open(v);
figure('Color','white')
ratio=5;%画晶格位置最后的放缩比例，使得画面不要超出figure
   hex_x = [0, 0.5*sqrt(3), 0.5*sqrt(3), 0, -0.5*sqrt(3), -0.5*sqrt(3),0] * a_cell;
    hex_y = [-1, -0.5, 0.5, 1, 0.5, -0.5,-1] * a_cell;
    hold on;axis equal;%set(gca, 'XColor', 'none', 'YColor', 'none');
    for i=1:Ny-1
        for j=1:(Nx-1)/2
        plot((hex_x+(j-1+(i-1)/2)*sqrt(3))/ratio, (hex_y-(i-1)*1.5)/ratio, 'LineWidth', 0.5,'Color',[0.8 0.8 0.8]); 
        end
    end
    U_new=0.01*ones(N_tot,nt);
    scatter_handle=scatter(save_num(:,2)/ratio, save_num(:,3)/ratio,abs(U_new(:,1)),angle(U_new(:,1)) , 'filled');patch([-1.1 NaN],[4 NaN],[-pi NaN],'EdgeColor','interp'),patch([-1.1 NaN],[4 NaN],[pi NaN],'EdgeColor','interp')
    colormap(cool);cbar = colorbar; title(cbar, 'phase');
    text_handle=text(0, 2, ['t=' num2str(t_vec(1))]);
    cbar.Ticks = [-pi,0, pi]; cbar.TickLabels = {'-\pi','0', '\pi'}; % Set the tick labels
for j=1:nt
    scatter_handle.SizeData = 30 * abs(phi(:,j))/max(abs(phi(:,j)));
            scatter_handle.CData = angle(phi(:,j)); 
            set(text_handle, 'String', ['t=' num2str(t_vec(j))]);
            axis off, pause(0.01)
            frame = getframe(gcf);writeVideo(v,frame)
end
    close(v);


for j=1:nt
    r_max2(j,1)=save_num(max_index(j),2)/sqrt(3);
    r_max2(j,2)=save_num(max_index(j),3)/sqrt(3);
end

figure('Color','white')
subplot(1,2,1), plot(k_max(:,1)/2/pi,k_max(:,2)/2/pi,'o','LineWidth',2),hold on
plot(k_max_the(:,1)/2/pi,k_max_the(:,2)/2/pi,'*','LineWidth',2)
title('k space'),xlabel('kx'),ylabel('ky'),ylim([0.4 0.5]),xlim([0.2 0.3])
%axis off
subplot(1,2,2),plot(r_max2(:,1),r_max2(:,2),'o','LineWidth',2),hold on
plot(r_max_the(:,1),r_max_the(:,2),'*','LineWidth',2)
title('r space'),xlabel('x'),ylabel('y')
%axis off
    legend('simulation','theory')
sgtitle('Maximum of Gaussian wave packet')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    L=301;
    k11=linspace(0,2*pi,L);k22=k11;w=[];
    for m=1:1:L
        for n=1:1:L
            k1=k11(m);k2=k22(n);
            w(m,n)=E_fun(k1,(k1+2*k2)/sqrt(3));
        end
    end
    [K1,K2]=meshgrid(k11,k22);
    KX=K1;KY=(K1+2*K2)/sqrt(3);


    xte=-0.3; yte=1;
figure('Color','white')
h1=subplot(4,2,1);pcolor(KX/2/pi,KY/2/pi,real(w)),cb=colorbar;
    xlabel('k_x'),ylabel('k_y'),title('Re(E)'),shading interp,axis equal
    xlim([0 1]),load('my_red.mat'),load('my_dark.mat');
    colormap(h1,my_red),cb.Position=[0.38    0.79    0.01    0.12];
    hold on,plot(KX(1,:)/2/pi,KY(1,:)/2/pi,'k--'),plot(KX(end,:)/2/pi,KY(end,:)/2/pi,'k--'),plot(KX(:,1)/2/pi,KY(:,1)/2/pi,'k--'),plot(KX(:,end)/2/pi,KY(:,end)/2/pi,'k--')
    text(xte-0.1, yte, '(a)', 'Units', 'normalized', 'FontSize', 12, 'HorizontalAlignment', 'center');
h2=subplot(4,2,2);pcolor(KX/2/pi,KY/2/pi,imag(w)),cb=colorbar;
    xlabel('k_x'),ylabel('k_y'),title('Im(E)'),shading interp,axis equal
    xlim([0 1]),colormap(h2,my_dark),cb.Position=[0.82    0.79    0.01    0.12];
    hold on,plot(KX(1,:)/2/pi,KY(1,:)/2/pi,'k--'),plot(KX(end,:)/2/pi,KY(end,:)/2/pi,'k--'),plot(KX(:,1)/2/pi,KY(:,1)/2/pi,'k--'),plot(KX(:,end)/2/pi,KY(:,end)/2/pi,'k--')
    text(xte-0.1, yte, '(b)', 'Units', 'normalized', 'FontSize', 12, 'HorizontalAlignment', 'center');
subplot('Position',[0.2 0.15,0.27 0.15])
plot(k_max(:,1)/2/pi,k_max(:,2)/2/pi,'o','LineWidth',2),hold on
    plot(k_max_the(:,1)/2/pi,k_max_the(:,2)/2/pi,'k','LineWidth',2)
    xlabel('k_x'),ylabel('k_y')%,ylim([0.42 0.49]),xlim([0.18 0.26])
    text(xte, yte, '(i)', 'Units', 'normalized', 'FontSize', 12, 'HorizontalAlignment', 'center');
    l=legend('simulation','theory');l.ItemTokenSize = [10,10];
subplot('Position',[0.63 0.15,0.27 0.15])
plot(r_max(:,1)-r_max(1,1),r_max(:,2)-r_max(1,2),'o','LineWidth',2),hold on
    plot(r_max_the(:,1)-r_max_the(1,1),r_max_the(:,2)-r_max_the(1,2),'r','LineWidth',2)
    xlabel('x'),ylabel('y'),xlim([0 1.5])
    text(xte+0.05, yte, '(j)', 'Units', 'normalized', 'FontSize', 12, 'HorizontalAlignment', 'center');
    l=legend('simulation','theory');l.ItemTokenSize = [10,10];

  h3=subplot(4,3,4);
     i=1;t=t_vec(i);
    pcolor(kx_grid/2/pi, ky_grid/2/pi, abs(X_K_T(:,:,i))); % 绘制x_k_t的绝对值
    shading interp; % 插值光滑
    cb=colorbar; axis equal,xlim([min(min(KX))/2/pi max(max(KX))/2/pi]),ylim([min(min(KY))/2/pi max(max(KY))/2/pi])
    xlabel('k_x'); ylabel('k_y'); % y轴标签
    text(0.2,1.6,['t = ' num2str(t)],'FontSize',12)
     hold on,plot(KX(1,:)/2/pi,KY(1,:)/2/pi,'k--'),plot(KX(end,:)/2/pi,KY(end,:)/2/pi,'k--'),plot(KX(:,1)/2/pi,KY(:,1)/2/pi,'k--'),plot(KX(:,end)/2/pi,KY(:,end)/2/pi,'k--')
    hold off,colormap(h3,my_sky_2),cb.Position=[0.324    0.57    0.01    0.12];
    text(xte, yte, '(c)', 'Units', 'normalized', 'FontSize', 12, 'HorizontalAlignment', 'center');

h4=subplot(4,3,5);
     i=9;t=t_vec(i);
    pcolor(kx_grid/2/pi, ky_grid/2/pi, abs(X_K_T(:,:,i))); % 绘制x_k_t的绝对值
    shading interp; % 插值光滑
    cb=colorbar; axis equal,xlim([min(min(KX))/2/pi max(max(KX))/2/pi]),ylim([min(min(KY))/2/pi max(max(KY))/2/pi])
    xlabel('k_x'); ylabel('k_y'); % y轴标签
    text(0.2,1.6,['t = ' num2str(t)],'FontSize',12)
     hold on,plot(KX(1,:)/2/pi,KY(1,:)/2/pi,'k--'),plot(KX(end,:)/2/pi,KY(end,:)/2/pi,'k--'),plot(KX(:,1)/2/pi,KY(:,1)/2/pi,'k--'),plot(KX(:,end)/2/pi,KY(:,end)/2/pi,'k--')
    hold off,colormap(h4,my_sky_2),cb.Position=[0.605    0.57    0.01    0.12];
    text(xte, yte, '(d)', 'Units', 'normalized', 'FontSize', 12, 'HorizontalAlignment', 'center');

h5=subplot(4,3,6);
     i=17;t=t_vec(i);
    pcolor(kx_grid/2/pi, ky_grid/2/pi, abs(X_K_T(:,:,i))); % 绘制x_k_t的绝对值
    shading interp; % 插值光滑
    cb=colorbar; axis equal,xlim([min(min(KX))/2/pi max(max(KX))/2/pi]),ylim([min(min(KY))/2/pi max(max(KY))/2/pi])
    xlabel('k_x'); ylabel('k_y'); % y轴标签
    text(0.2,1.6,['t = ' num2str(t)],'FontSize',12)
     hold on,plot(KX(1,:)/2/pi,KY(1,:)/2/pi,'k--'),plot(KX(end,:)/2/pi,KY(end,:)/2/pi,'k--'),plot(KX(:,1)/2/pi,KY(:,1)/2/pi,'k--'),plot(KX(:,end)/2/pi,KY(:,end)/2/pi,'k--')
    hold off,colormap(h5,my_sky_2),cb.Position=[0.885    0.57    0.01    0.12];
    text(xte, yte, '(e)', 'Units', 'normalized', 'FontSize', 12, 'HorizontalAlignment', 'center');



subplot(4,3,7)
ratio=5;%画晶格位置最后的放缩比例，使得画面不要超出figure
   hex_x = [0, 0.5*sqrt(3), 0.5*sqrt(3), 0, -0.5*sqrt(3), -0.5*sqrt(3),0] * a_cell;
    hex_y = [-1, -0.5, 0.5, 1, 0.5, -0.5,-1] * a_cell;
    hold on;axis equal;%set(gca, 'XColor', 'none', 'YColor', 'none');
    for i=1:Ny-1
        for j=1:(Nx-1)/2
        plot((hex_x+(j-1+(i-1)/2)*sqrt(3))/ratio, (hex_y-(i-1)*1.5)/ratio, 'LineWidth', 0.5,'Color',[0.8 0.8 0.8]); 
        end
    end
    U_new=0.01*ones(N_tot,nt);
    scatter_handle=scatter(save_num(:,2)/ratio, save_num(:,3)/ratio,abs(U_new(:,1)),'b' , 'filled');patch([-1.1 NaN],[4 NaN],[-pi NaN],'EdgeColor','interp'),patch([-1.1 NaN],[4 NaN],[pi NaN],'EdgeColor','interp')
    j=1;
    scatter_handle.SizeData = 15 * abs(phi(:,j))/max(abs(phi(:,j)));
    text(6, -3, ['t=' num2str(t_vec(j))],'FontSize',12);       
    axis off
           xlim([5.3,11.6]),ylim([-8.2,-2])
           text(xte+0.3, yte, '(f)', 'Units', 'normalized', 'FontSize', 12, 'HorizontalAlignment', 'center');

subplot(4,3,8)
ratio=5;%画晶格位置最后的放缩比例，使得画面不要超出figure
   hex_x = [0, 0.5*sqrt(3), 0.5*sqrt(3), 0, -0.5*sqrt(3), -0.5*sqrt(3),0] * a_cell;
    hex_y = [-1, -0.5, 0.5, 1, 0.5, -0.5,-1] * a_cell;
    hold on;axis equal;%set(gca, 'XColor', 'none', 'YColor', 'none');
    for i=1:Ny-1
        for j=1:(Nx-1)/2
        plot((hex_x+(j-1+(i-1)/2)*sqrt(3))/ratio, (hex_y-(i-1)*1.5)/ratio, 'LineWidth', 0.5,'Color',[0.8 0.8 0.8]); 
        end
    end
    U_new=0.01*ones(N_tot,nt);
    scatter_handle=scatter(save_num(:,2)/ratio, save_num(:,3)/ratio,abs(U_new(:,1)),'b' , 'filled');patch([-1.1 NaN],[4 NaN],[-pi NaN],'EdgeColor','interp'),patch([-1.1 NaN],[4 NaN],[pi NaN],'EdgeColor','interp')
    j=9;
    scatter_handle.SizeData = 15 * abs(phi(:,j))/max(abs(phi(:,j)));
    text(6, -3, ['t=' num2str(t_vec(j))],'FontSize',12);     
    axis off
           xlim([5.3,11.6]),ylim([-8.2,-2])
           text(xte+0.3, yte, '(g)', 'Units', 'normalized', 'FontSize', 12, 'HorizontalAlignment', 'center');



subplot(4,3,9)
ratio=5;%画晶格位置最后的放缩比例，使得画面不要超出figure
   hex_x = [0, 0.5*sqrt(3), 0.5*sqrt(3), 0, -0.5*sqrt(3), -0.5*sqrt(3),0] * a_cell;
    hex_y = [-1, -0.5, 0.5, 1, 0.5, -0.5,-1] * a_cell;
    hold on;axis equal;%set(gca, 'XColor', 'none', 'YColor', 'none');
    for i=1:Ny-1
        for j=1:(Nx-1)/2
        plot((hex_x+(j-1+(i-1)/2)*sqrt(3))/ratio, (hex_y-(i-1)*1.5)/ratio, 'LineWidth', 0.5,'Color',[0.8 0.8 0.8]); 
        end
    end
   j=17;
    scatter(save_num(:,2)/ratio, save_num(:,3)/ratio,15 * abs(phi(:,j))/max(abs(phi(:,j))),'b' , 'filled');patch([-1.1 NaN],[4 NaN],[-pi NaN],'EdgeColor','interp'),patch([-1.1 NaN],[4 NaN],[pi NaN],'EdgeColor','interp')
   text(6, -3, ['t=' num2str(t_vec(j))],'FontSize',12);
            axis off
           xlim([5.3,11.6]),ylim([-8.2,-2])
           text(xte+0.3, yte, '(h)', 'Units', 'normalized', 'FontSize', 12, 'HorizontalAlignment', 'center');
