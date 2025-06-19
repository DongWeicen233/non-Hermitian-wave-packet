function [save_num,marker,y_unique] = coordinate(Nx,Ny)
% Author: Weicen Dong <orange233@sjtu.edu.cn>
% Created Date: 2025/1/21

    %save_num 保存原子的编号、x、y坐标
    %marker 保存原子的标记值,第j个元素值为i表示第j个原子的y坐标值为y_unique(i)，从大到小排第i个
    %y_unique 保存所有不同的y坐标值
    pl = 1;
    for ny = 1:Ny
        for nx = 1:Nx+1
            if (ny == 1 || ny == Ny) && nx == Nx+1
                break;
            else
                [N, x, y] = position(Nx, ny, nx);
                save_num(pl,1) = N;save_num(pl,2) = x;save_num(pl,3) = y;
                pl = pl+ 1;
            end
        end
    end

% 按 y 坐标从大到小排序 save_num 矩阵
save_num = sortrows(save_num, -3);  % 第3列是 y 坐标，按 y 从大到小排序
y_unique = flip(unique(save_num(:,3)));  % 获取所有不同的 y 坐标值
marker = zeros(size(save_num,1),1);  % 初始化标记数组

for i = 1:length(y_unique)
    y_val = y_unique(i);  % 当前 y 坐标值
    % 查找具有相同 y 值的原子
    idx = find(save_num(:,3) == y_val);
    marker(idx) = i;  % 对应原子的标记值为 i
end

end