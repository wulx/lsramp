function varargout = lsr_plot(f_list, dt_list)
%LSR_PLOT Linear Speed Ramp PLOT
%
% varargin:
%   f_list   --  frequencies list
%   dt_list  --  time periods list
%   t_step   --  time step
%
% varargout:
%   p_a   --   polynomial curve fitting outputs for acceleration ramp
%   S_a   --   ditto
%   p_d   --   for deceleration ramp
%   S_d   --   ditto

% copyright (c) wulx, <gurdy.woo@mail.ustc.edu.cn>
% last modified by wulx, 2013/10/31

sn = numel(dt_list);
t_list = [0 arrayfun(@(n) sum(dt_list(1:n)), 1:sn)];
p_list = [f_list f_list(end)];

figure, hold on;
stairs(t_list, p_list, 'k-')
xlim([0 t_list(end)])

% baseline
plot([0 t_list(end)], [0 0], 'k-')

% find the max frequencies
idx = find(f_list == max(f_list));

% speeding up ramp
xc_up = arrayfun(@(i) t_list(i)+0.5*dt_list(i), 1:idx(1));
yc_up = f_list(1:idx(1));
plot(xc_up, yc_up, 'r.')

[p_a, S_a] = polyfit(xc_up, yc_up, 1); % Polynomial curve fitting for acceleration ramp
xc_poly = [0 xc_up t_list(idx(2))];
yc_poly = polyval(p_a, xc_poly);

plot(xc_poly, yc_poly, 'k:')
plot(xc_poly([end end]), [0 yc_poly(end)], 'k:')

% speeding down ramp
xc_dn = arrayfun(@(i) t_list(i)+0.5*dt_list(i), idx(end):sn);
yc_dn = f_list(idx(end):sn);
plot(xc_dn, yc_dn, 'r.')

[p_d, S_d] = polyfit(xc_dn, yc_dn, 1); % for deceleration ramp
xc_poly2 = [t_list(idx(end)) xc_dn t_list(end)];
yc_poly2 = polyval(p_d, xc_poly2);

plot(xc_poly2, yc_poly2, 'k:')
plot(xc_poly2([1 1]), [0 yc_poly2(1)], 'k:')

% set Y axis limit adaptely
ylim([min([yc_poly(1), yc_poly2(end), 0]), max(yc_poly(end), yc_poly2(1))])

title(['Acceleration: ' num2str(p_a(1)) ' ; deceleration: ' num2str(p_d(1))]);

switch nargout
    case 2
        varargout = {p_a, p_d};
    case 4
        varargout = {p_a, S_a, p_d, S_d};
    otherwise
        disp('well done:)');
end
