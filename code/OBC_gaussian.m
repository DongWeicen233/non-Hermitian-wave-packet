% Author: Weicen Dong <orange233@sjtu.edu.cn>
% Created Date: 2025/1/21
% 
% %高斯波包在strip上移动
clear all

load('data_OBC_NH.mat')
k_vec=kk(635:900)/2/pi;e_vec=ek(n,635:900);
for i=635:900,U_vec(:,i-634)=U(:,n,i);end


% 使用样条插值,之后需要用到能量求导得到理论上的x_max和k——max
k_fine = linspace(min(k_vec), max(k_vec), 1000); % 在 k_vec 范围内生成更细的 k 点
e_interp = spline(k_vec, e_vec, k_fine); % 插值得到光滑的 e(k)
e_k_fun = @(k) interp1(k_fine,e_interp,k, 'linear', 'extrap');
de_dk = gradient(e_interp, k_fine); % 计算插值后的 e(k) 的导数
de_dk_func = @(k) interp1(k_fine, de_dk, k, 'linear', 'extrap');

%gaussian wave packet in k space
%初始化波包，时间矢量，等空矩阵
k_bar=0.45;Delta=0.02;
%F=-0.0005;
F=0;

x_k_0=@(k,t) exp(-((k-k_bar-F*t)/Delta).^2/2)/sqrt(2*pi)/Delta;
psi_k(1,:)=x_k_0(k_vec,0);
nt=40;t_vec=linspace(0,20,nt);dt=t_vec(2)-t_vec(1);

k_max_the=zeros(nt,1);x_max_the=zeros(nt,1);
k_max_the(1)=k_bar;x_max_the(1)=30;


%高斯波包在k空间的演化
v = VideoWriter(sprintf('video_gaussian_k.mp4'), 'MPEG-4');
rate=5; v.FrameRate = rate;open(v);
Et=@(k,t) arrayfun(@(k) trapz(linspace(0,t,ceil(t*3)),e_k_fun(k+F*linspace(0,t,ceil(t*3)))),k);

figure('Color','w')
for i=1:nt
    t=t_vec(i);
    psi_k(i,:)=x_k_0(k_vec,t).*exp(-1i*Et(k_vec,t));%波包随时间演化，不同k的波包随时间演化的相位不同
    sgtitle(['t=' num2str(t_vec(i))])
    plot(k_vec,abs(psi_k(i,:)),'o'),ylim([0 1.2*max(abs(psi_k(i,:)))]),xlabel('k'),ylabel('|\psi_k|')
    frame = getframe(gcf);writeVideo(v,frame);
    if i>1
    k_max_the(i)=F*dt+k_max_the(i-1)+Delta^2*dt*imag(de_dk_func(k_max_the(i-1)));
    x_max_the(i)=x_max_the(i-1)+dt*real(de_dk_func(k_max_the(i-1)))/2/pi;%多除了2pi是因为k的周期实际上是2pi。算k时没除是因为我们把k的周期当成1
    end
end
close(v);

%通过傅立叶变换得到实空间的分布
Ny=n;Nx=151;N_tot=Nx*2+(Ny-2)*(Nx+1);a_cell=1;[save_num,marker,y_unique] = coordinate(Nx,Ny);
v = VideoWriter(sprintf('video_gaussian_r.mp4'), 'MPEG-4');
rate=5; v.FrameRate = rate;open(v);

figure('Color', 'w'),ratio=5;%画晶格位置最后的放缩比例，使得画面不要超出figure
   hex_x = [0, 0.5*sqrt(3), 0.5*sqrt(3), 0, -0.5*sqrt(3), -0.5*sqrt(3),0] * a_cell;
    hex_y = [-1, -0.5, 0.5, 1, 0.5, -0.5,-1] * a_cell;
    hold on;axis equal;set(gca, 'XColor', 'none', 'YColor', 'none');
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
for i=1:nt%不同的时刻用的psi_k不同
    for pp=1:N_tot%计算每个原子的U
        U_new(pp,i)=sum(psi_k(i,:).*U_vec(marker(pp),:).*exp(1i*(save_num(pp,2)/sqrt(3)-x_max_the(1))*k_vec*2*pi));%根据第pp个点的x坐标算相位,k_vec是除了2pi的，现在要乘上去
        %sum是对每个k求和
        %marker(pp)表示第pp个原子的y坐标是从上到下多少个
        %这里用的x-20就是说波包的中心在20，x/sqrt(3)之后是以我的a1为基矢的横坐标
    end
    scatter_handle.SizeData = 30 * abs(U_new(:, i))/max(abs(U_new(:, i)));
        scatter_handle.CData = angle(U_new(:, i)); 
        set(text_handle, 'String', ['t=' num2str(t_vec(i))]);
        axis off,pause(0.01)
        frame = getframe(gcf);writeVideo(v,frame);
end
close(v);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%画出计算和理论的xk最大值以及能带
%对每个t求得U_new最大的位置,然后得到这个位置对应的x，再画出x-t图
for i=1:nt,[max_val,max_index]=max(abs(U_new(:,i)));x_max(i)=save_num(max_index,2);end
for i=1:nt,[max_val,max_index]=max(abs(psi_k(i,:)));k_max(i)=k_vec(max_index);end
figure('Color','w')
subplot(2,2,1),plot(t_vec,x_max/sqrt(3),'o'),xlabel('t'),ylabel('x_max'),hold on
plot(t_vec,x_max_the,'LineWidth',2),legend('simulation','theory')
subplot(2,2,2),plot(t_vec,k_max,'o'),xlabel('t'),ylabel('k_max'),hold on
plot(t_vec,k_max_the,'LineWidth',2),legend('simulation','theory')
subplot(2,2,3),plot(k_vec,real(e_vec)),xlabel('k'),ylabel('Re(e_k)')
subplot(2,2,4),plot(k_vec,imag(e_vec)),xlabel('k'),ylabel('Im(e_k)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Color','w'),xte=-0.1; yte=1.15;
subplot(4,2,8),plot(t_vec,x_max/sqrt(3)-30,'o'),xlabel('t'),ylabel('x_{M}'),hold on
plot(t_vec,x_max_the-30,'LineWidth',2,'Color','r'),l=legend('simulation','theory');l.ItemTokenSize = [10,10];
text(xte, yte, '(e)', 'Units', 'normalized', 'FontSize', 12, 'HorizontalAlignment', 'center');
subplot(4,2,7),plot(t_vec,k_max,'o'),xlabel('t'),ylabel('k_{M}'),hold on
plot(t_vec,k_max_the,'LineWidth',2,'Color','k'),l=legend('simulation','theory');l.ItemTokenSize = [10,10];
text(xte, yte, '(d)', 'Units', 'normalized', 'FontSize', 12, 'HorizontalAlignment', 'center');

%subplot(4,2,3,'Position',[0.15 0.55 0.7 0.15]);
subplot(4,2,2)
yyaxis left; 
plot(k_vec,abs(psi_k(1,:)),'b','LineWidth',2),ylim([0 1.2*max(abs(psi_k(1,:)))]),xlabel('k'),ylabel('|\psi_k|'),set(gca, 'YColor', 'b');
yyaxis right;
plot(k_vec,abs(psi_k(end,:)),'m','LineWidth',2),ylim([0 1.2*max(abs(psi_k(end,:)))]),xlabel('k'),ylabel('|\psi_k|'),set(gca, 'YColor', 'm');
box on,xlim([0.35 0.5])
ee1=e_k_fun(k_max_the(1));ee20=e_k_fun(k_max_the(end));
l=legend('t=0','t=20');l.ItemTokenSize = [10,10];
text(xte, yte, '(b)', 'Units', 'normalized', 'FontSize', 12, 'HorizontalAlignment', 'center');
% text(0.82,0.5,['E_e=' num2str(ee1,'%0.2f')], 'Units', 'normalized', 'FontSize', 10, 'HorizontalAlignment', 'center','Color','b')
% text(0.18,0.82,['E_e=' num2str(ee20,'%0.2f')], 'Units', 'normalized', 'FontSize', 10, 'HorizontalAlignment', 'center','Color','m')

subplot(4,2,1)
yyaxis left; 
plot(k_vec,real(e_vec),'r','LineWidth',2),ylabel('Re(E_k)'),set(gca, 'YColor', 'r');
yyaxis right;
plot(k_vec,imag(e_vec),'k','LineWidth',2),xlabel('k'),ylabel('Im(E_k)'),set(gca, 'YColor', 'k');
box on
text(xte, yte, '(a)', 'Units', 'normalized', 'FontSize', 12, 'HorizontalAlignment', 'center');


subplot(4,1,2)
hex_x = [0, 0.5*sqrt(3), 0.5*sqrt(3), 0, -0.5*sqrt(3), -0.5*sqrt(3),0] * a_cell;
hex_y = [-1, -0.5, 0.5, 1, 0.5, -0.5,-1] * a_cell;
hold on;axis equal;set(gca, 'XColor', 'none', 'YColor', 'none');
for i=1:Ny-1
    for j=1:(Nx-1)/2
    plot((hex_x+(j-1+(i-1)/2)*sqrt(3))/ratio, (hex_y-(i-1)*1.5)/ratio, 'LineWidth', 0.5,'Color',[0.8 0.8 0.8]); 
    end
end
U_new=0.01*ones(N_tot,nt);
scatter_handle=scatter(save_num(:,2)/ratio, save_num(:,3)/ratio,abs(U_new(:,1)),'b' , 'filled');patch([-1.1 NaN],[4 NaN],[-pi NaN],'EdgeColor','interp'),patch([-1.1 NaN],[4 NaN],[pi NaN],'EdgeColor','interp')
colormap(cool);
i=1;
for pp=1:N_tot%计算每个原子的U
    U_new(pp,i)=sum(psi_k(i,:).*U_vec(marker(pp),:).*exp(1i*(save_num(pp,2)/sqrt(3)-x_max_the(1))*k_vec*2*pi));%根据第pp个点的x坐标算相位,k_vec是除了2pi的，现在要乘上去
    %sum是对每个k求和
    %marker(pp)表示第pp个原子的y坐标是从上到下多少个
    %这里用的x-20就是说波包的中心在20，x/sqrt(3)之后是以我的a1为基矢的横坐标
end
scatter_handle.SizeData = 20 * abs(U_new(:, i))/max(abs(U_new(:, i)));
    %set(text_handle, 'String', ['t=' num2str(t_vec(i))]);
    axis off,pause(0.01)
    % cbar = colorbar; title(cbar, 'phase');
    % text_handle=text(0, 2, ['t=' num2str(t_vec(1))]);
    % cbar.Ticks = [-pi,0, pi]; cbar.TickLabels = {'-\pi','0', '\pi'}; % Set the tick labels
    xlim([4 20.1]),ylim([-5.6 0.2])
    text(0.1, 0.9, 't=0', 'Units', 'normalized', 'FontSize', 12, 'HorizontalAlignment', 'center');
    text(xte+0.1, yte-0.2, '(c)', 'Units', 'normalized', 'FontSize', 12, 'HorizontalAlignment', 'center');

    subplot(4,1,3)
    hex_x = [0, 0.5*sqrt(3), 0.5*sqrt(3), 0, -0.5*sqrt(3), -0.5*sqrt(3),0] * a_cell;
    hex_y = [-1, -0.5, 0.5, 1, 0.5, -0.5,-1] * a_cell;
    hold on;axis equal;set(gca, 'XColor', 'none', 'YColor', 'none');
    for i=1:Ny-1
        for j=1:(Nx-1)/2
        plot((hex_x+(j-1+(i-1)/2)*sqrt(3))/ratio, (hex_y-(i-1)*1.5)/ratio, 'LineWidth', 0.5,'Color',[0.8 0.8 0.8]); 
        end
    end
    U_new=0.01*ones(N_tot,nt);
    scatter_handle=scatter(save_num(:,2)/ratio, save_num(:,3)/ratio,abs(U_new(:,1)),'m' , 'filled');patch([-1.1 NaN],[4 NaN],[-pi NaN],'EdgeColor','interp'),patch([-1.1 NaN],[4 NaN],[pi NaN],'EdgeColor','interp')
    colormap(cool);
    i=nt;
    for pp=1:N_tot%计算每个原子的U
        U_new(pp,i)=sum(psi_k(i,:).*U_vec(marker(pp),:).*exp(1i*(save_num(pp,2)/sqrt(3)-x_max_the(1))*k_vec*2*pi));%根据第pp个点的x坐标算相位,k_vec是除了2pi的，现在要乘上去
        %sum是对每个k求和
        %marker(pp)表示第pp个原子的y坐标是从上到下多少个
        %这里用的x-20就是说波包的中心在20，x/sqrt(3)之后是以我的a1为基矢的横坐标
    end
    scatter_handle.SizeData = 20 * abs(U_new(:, i))/max(abs(U_new(:, i)));

        %set(text_handle, 'String', ['t=' num2str(t_vec(i))]);
        axis off,pause(0.01)
        % cbar = colorbar; title(cbar, 'phase');
        % text_handle=text(0, 2, ['t=' num2str(t_vec(1))]);
        % cbar.Ticks = [-pi,0, pi]; cbar.TickLabels = {'-\pi','0', '\pi'}; % Set the tick labels
        xlim([4 20.1]),ylim([-5.6 0.2])
        text(0.1, 0.9, 't=20', 'Units', 'normalized', 'FontSize', 12, 'HorizontalAlignment', 'center');



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %实空间哈密顿量作用在波包上
% U_evol(:,1)=zeros(length(U_new(:,1)),1); U_evol(10)=0.1;


% [E,V,H] = Hamiltonian(Nx,Ny,t1,a,b,c,d);
% nt=20;t_vec=linspace(0,20,nt);dt=t_vec(2)-t_vec(1);

% figure('Color', 'w'),ratio=5;
%    hex_x = [0, 0.5*sqrt(3), 0.5*sqrt(3), 0, -0.5*sqrt(3), -0.5*sqrt(3),0] * a_cell;
%     hex_y = [-1, -0.5, 0.5, 1, 0.5, -0.5,-1] * a_cell;
%     hold on;axis equal;set(gca, 'XColor', 'none', 'YColor', 'none');
%     for i=1:Ny-1
%         for j=1:(Nx-1)/2
%         plot((hex_x+(j-1+(i-1)/2)*sqrt(3))/ratio, (hex_y-(i-1)*1.5)/ratio, 'LineWidth', 0.5,'Color',[0.8 0.8 0.8]); 
%         end
%     end
% %scatter_handle=scatter(save_num(:,2)/ratio, save_num(:,3)/ratio,1 * abs(U_new(:, 1))/max(abs(U_new(:, 1))),angle(U_new(:,1)) , 'filled');patch([-1.1 NaN],[4 NaN],[-pi NaN],'EdgeColor','interp'),patch([-1.1 NaN],[4 NaN],[pi NaN],'EdgeColor','interp')
% scatter_handle=scatter(save_num(:,2)/ratio, save_num(:,3)/ratio,1 * abs(U_new(:, 1))/max(abs(U_new(:, 1))),'b' , 'filled');patch([-1.1 NaN],[4 NaN],[-pi NaN],'EdgeColor','interp'),patch([-1.1 NaN],[4 NaN],[pi NaN],'EdgeColor','interp')

% colormap(cool);cbar = colorbar; title(cbar, 'phase');
% text_handle=text(0, 2, ['t=' num2str(t_vec(1))]);
% cbar.Ticks = [-pi,0, pi]; cbar.TickLabels = {'-\pi','0', '\pi'}; % Set the tick labels

% v = VideoWriter(sprintf('video_gaussian_square.mp4'), 'MPEG-4');
% rate=5; v.FrameRate = rate;open(v);
% for i=1:nt
% U_evol(:,i)=expm(-1i*H*t_vec(i))*U_new(:,1);
% scatter_handle.SizeData = 30 * abs(U_evol(:, i))/max(abs(U_evol(:, i)));
% scatter_handle.CData = angle(U_evol(:, i)); 
% set(text_handle, 'String', ['t=' num2str(t_vec(i))]);
% axis off,pause(0.01)
% frame = getframe(gcf);writeVideo(v,frame);
% end
% close(v)
