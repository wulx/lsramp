function [f_list, dt_list] = steps_per_time(sn_a, sn_c, sn_d, pf_i, pf_m, s_u, method)
%STEPS_PER_TIME steps per time algorithm
%
% varargin:
%   sn_a   --   stepper numbers of acceleration
%   sn_c   --   stepper numbers during moving with constant speed
%   sn_d   --   stepper numbers of deceleration
%   pf_i   --   initial pulse frequency
%   pf_m   --   maximum pulse frequency
%   s_u    --   unit steps, default value is 1
%   method --   approximation method
%
% varargout:
%   f_list   --  frequencies list
%   dt_list  --  time periods list

% copyright (c) wulx, <gurdy.woo@gmail.com>
% last modified by wulx, 2013/10/16

% default settings
if nargin < 7, method = 'ideal'; end
if nargin < 6, s_u = 1; end

sn = sn_a + sn_d - 1;
f_list = zeros(1, sn); % pre-allocation for better performance
dt_list = zeros(1, sn); % ditto

% unit time period
dt = s_u / pf_i;

switch method
    case 'ideal'
        dt_list(1:end) = dt;
        % #1 acceleration
        % divide frequencies uniformly in accordance with frequency range (pf_i, pf_m) and step numbers (sn_a)
        f_list(1:sn_a) = linspace(pf_i, pf_m, sn_a);
        
        % #2 moving with constant speed
        dt_list(sn_a) = (sn_c + 2) * dt;
        
        % #3 deceleration
        % set f_list(sn_a) to pf_m again, duplicated but harmless!
        f_list(sn_a - 1 + (1:sn_d)) = linspace(pf_m, pf_i, sn_d); 
    case {'round', 'fix'}
        % select round method
        %round_to = @round; % function handle
        %if strcmp(method, 'fix')
        %    round_to = @fix;
        %end
        round_to = @(x) feval(method, x); % more concise
        
        % #1 acceleration
        % divide frequencies uniformly and round them to nearest integer
        f_list1 = round_to( linspace(pf_i, pf_m, sn_a) );
        s_list1 = round_to( dt * f_list1); % step or pulse numbers should be integer
        dt_list1 = s_list1 ./ f_list1; % actual time periods
        
        f_list(1:sn_a-1) = f_list1(1:end-1); % from 1 to sn_a-1
        dt_list(1:sn_a-1) = dt_list1(1:end-1);
        
        % #2 moving with constant speed
        s_tot2 = round_to( (sn_c + 2) * pf_m * dt ); % total step numbers (rounded to interger)
        dt2 = s_tot2 / pf_m;
        
        f_list(sn_a) = pf_m; % sn_a
        dt_list(sn_a) = dt2;
        
        % #3 deceleration
        f_list3 = round_to( linspace(pf_m, pf_i, sn_d) );
        s_list3 = round_to( dt * f_list3 );
        dt_list3 = s_list3 ./ f_list3;
        
        inds = sn_a + (1:sn_d-1); % from sn_a+1 to sn_a+sn_d-1 (i.e. sn)
        f_list(inds) = f_list3(2:end);
        dt_list(inds) = dt_list3(2:end);
    otherwise
        disp(['unknown approximation method: ' method])
        disp('available methods: ideal(default), round and fix')
end
