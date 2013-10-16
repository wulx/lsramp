function lsr_plot(f_list, dt_list, t_step)
%LSR_PLOT Linear Speed Ramp PLOT
%
% varargin:
%   f_list   --  frequencies list
%   dt_list  --  time periods list
%   t_step   --  time step
%
% varargout:
%   t_seq   --  timeline sequences
%   f_seq   --  sequence of frequencies
%   s_seq   --  sequence of accumulative stepper numbers
%
% copyright (c) wulx, <gurdy.woo@mail.ustc.edu.cn>
% last modified by wulx, 2013/10/16

if nargin < 3
    t_step = min( dt_list ) / 5;
end

t_tot = sum( dt_list ); % total elapsed time
t_seq = linspace(0, t_tot, round(t_tot / t_step));

f_seq = zeros( size(t_seq) );

num = numel(dt_list);
t_pts = arrayfun(@(n) sum( dt_list(1:n) ), 0:num); % time points

% piecewise function
for i = 1:num-1
    f_seq = f_list(i) * xor(t_seq < t_pts(i), t_seq < t_pts(i+1)) + f_seq;
end
f_seq = f_list(num) * xor(t_seq < t_pts(num), t_seq <= t_pts(num+1)) + f_seq;

figure, hold on;
plot(t_seq, f_seq, 'b-')

center_pts = arrayfun(@(n) t_pts(n) + 0.5*dt_list(n), 1:num); % all center time points
plot(center_pts, f_list, 'r.')

[~, ind] = max(f_list); % max frequency as watershed
acc_pts = center_pts(1:ind-1); % acceleration time points
dec_pts = center_pts((ind+1):end); % deceleration time points
f_acc = f_list(1:ind-1);
f_dec = f_list((ind+1):end);

plot(acc_pts, f_acc, dec_pts, f_dec, 'k-')

p_a = polyfit(acc_pts, f_acc, 1); % Polynomial curve fitting for acceleration ramp
p_d = polyfit(dec_pts, f_dec, 1); % for deceleration ramp

title(['Acceleration: ' num2str(p_a(1)) ' ; deceleration: ' num2str(p_d(1))]);

