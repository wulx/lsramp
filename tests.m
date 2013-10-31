%% testing on lsramp
clear all; close all; clc

%% time_per_step
sn_a = 30;
sn_c = 10;
sn_d = 18;
pf_i = 100;
pf_m = 1000;
s_u = 1;
method = 'round';
[f_list, dt_list] = time_per_step(sn_a, sn_c, sn_d, pf_i, pf_m, s_u, method);
% 
% sn = numel(dt_list);
% t_list = [0 arrayfun(@(n) sum(dt_list(1:n)), 1:sn)];
% f_list2 = [f_list f_list(end)];
% figure, hold on;
% stairs(t_list, f_list2)
% xlim([0 t_list(end)])
% ylim([0 pf_m])
% 
% [xb, yb] = stairs(t_list, f_list2);
% 
% idx = find(f_list == max(f_list));
% % speed up
% xUp = arrayfun(@(i) t_list(i)+0.5*dt_list(i), 1:idx(1));
% plot(xUp, f_list(1:idx(1)), 'r.')

lsr_plot(f_list, dt_list)
sn_plot(f_list, dt_list)

%% steps_per_time
% [f_list, dt_list] = steps_per_time(sn_a, sn_c, sn_d, pf_i, pf_m, s_u, method);
% 
% sn = numel(dt_list);
% t_list = [0 arrayfun(@(n) sum(dt_list(1:n)), 1:sn)];
% f_list2 = [f_list f_list(end)];
% figure, hold on;
% stairs(t_list, f_list2)
% xlim([0 t_list(end)])
% ylim([0 pf_m])
% 
% stairs(t_list, f_list2, 'ro')
% 
% idx = find(f_list == max(f_list));
% % speed up
% xUp = arrayfun(@(i) t_list(i)+0.5*dt_list(i), 1:sn);
% plot(xUp, f_list(1:sn), 'r.')
