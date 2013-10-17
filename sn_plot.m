function [s_seq, t_seq] = sn_plot(f_list, dt_list)
%SN_PLOT Stepper Numbers (accumulative) PLOT
%
% varargin:
%   f_list   --  frequencies list
%   dt_list  --  time periods list
%
% varargout:
%
% copyright (c) wulx, <gurdy.woo@gmail.com>
% last modified by wulx, 2013/10/17

num = numel(f_list); % number of frequencies
% fix bug #1 add round
sn_list = round( f_list .* dt_list ); % stepper numbers (of every frequency) list
sn_tot = sum( sn_list ); % total stepper numbers

t_list = arrayfun(@(n) sum( dt_list(1:n) ), 0:num);
t_seq = zeros(1, sn_tot+1);

for i = 1:num
    sn_i = sn_list(i);
    sn_a = sum( sn_list(1:i-1) );
    t_seq(sn_a + (1:sn_i)) = t_list(i) + (0:sn_i-1) / f_list(i);
end
t_seq(end) = t_seq(end-1); % only sn_tot pulses

s_seq = [1:sn_tot sn_tot];

figure, hold on;
plot(t_seq, s_seq, 'b-')
plot(t_seq, s_seq, 'r.')
