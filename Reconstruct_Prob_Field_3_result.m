clc, clear, close all;
transition_rates_xr = importdata('transition_rates_x+.mat');
transition_rates_xl = importdata('transition_rates_x-.mat');
transition_rates_yr = importdata('transition_rates_y+.mat');
transition_rates_yl = importdata('transition_rates_y-.mat');
%% ------将4个方向transition_rates_xr重组，得到重构的向量平均场------------
%%
% 初始化，注意与 "Reconstruct_Prob_Field_2_main" 保持一致
grid_size = [4, 4];  % 设置网格大小
bundary_x_l = 100; bundary_x_r = 200;  % 设置 x 轴边界（统计全局）
bundary_y_l = 100; bundary_y_r = 200;  % 设置 y 轴边界（统计全局） 
num_grids = [(bundary_x_r - bundary_x_l)/grid_size(1), (bundary_y_r - bundary_y_l)/grid_size(2)];  

%%
% 重组,整合,作图
transition_rates_x = transition_rates_xr - transition_rates_xl;
transition_rates_x = transition_rates_x(2:end-1, 2:end-1);
transition_rates_y = transition_rates_yr - transition_rates_yl;
transition_rates_y = transition_rates_y(2:end-1, 2:end-1);

% 计算网格中心点的坐标
[grid_x, grid_y] = meshgrid(bundary_x_l + grid_size(1)/2 : grid_size(1) : bundary_x_r, bundary_y_l + grid_size(2)/2 : grid_size(2) : bundary_y_r);
grid_x = grid_x(2:end-1, 2:end-1)';
grid_y = grid_y(2:end-1, 2:end-1)';

% 绘制向量场
% quiver(grid_x, grid_y, transition_rates_x, zeros(size(transition_rates_x)));
% quiver(grid_x, grid_y, zeros(size(transition_rates_y)), transition_rates_y);
h = quiver(grid_x, grid_y, transition_rates_x, transition_rates_y, 1.5);
set(h, 'LineWidth', 1);   % 调整向量场的粗细
xticks(bundary_x_l:10:bundary_x_r);   % 设置 x 轴和 y 轴的刻度
yticks(bundary_y_l:10:bundary_y_r);
xlabel('x_1');
ylabel('x_2');
title('Transition Prob. Rate (Mean) Field');