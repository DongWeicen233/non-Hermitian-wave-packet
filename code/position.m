function [N,x,y] = position(Nx,ny,nx)
% Author: Weicen Dong <orange233@sjtu.edu.cn>
% Created Date: 2025/1/21

%input
% Nx Ny 横纵有多少点
%nx ny 某个点的位置，｜ny,nx>
%output N, 是第几个点
if ny>1
N=nx+(ny-2)*(Nx+1)+Nx;
else 
    N=nx;
end 

%六边形边长为1，左上角六边形中心是原点
if ny==1
x=(nx-2)*sqrt(3)/2;
else 
x=(nx-2+(ny-2))*sqrt(3)/2;
end

y=1-(ny-1)*1.5-mod(N,2)*0.5;

end
