function varargout = sn_plot(f_list, dt_list)
%SN_PLOT number of steps vs time PLOT
%
% varargin:
%   f_list   --  frequencies list
%   dt_list  --  time periods list
%
% varargout:
%   s_seq   --   accumulative steps sequence
%   t_seq   --   timeline

% copyright (c) wulx, <gurdy.woo@gmail.com>
% last modified by wulx, 2013/10/31

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
t_seq(end) = t_seq(end-1) + 1/f_list(end); % no pulse

s_seq = [1:sn_tot sn_tot];

switch nargout
    case 0
        % plot stepping profile when no output arguments
        figure, hold on;
        plot(t_seq, s_seq, 'k-')
        plot(t_seq, s_seq, 'r.')
        
        xlim(t_list([1 end]))
    case 1
        varargout = {s_seq};
    case 2
        varargout = {s_seq, t_seq};
    otherwise
        error('number of output arguments should be less than 3.')
end
