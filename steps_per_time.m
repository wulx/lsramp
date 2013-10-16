function [f_list, dt_list] = steps_per_time(sn_a, sn_c, sn_d, pf_i, pf_m, s_u)
%STEPS_PER_TIME steps per time algorithm
%
% varargin:
%   sn_a   --   stepper numbers of acceleration
%   sn_c   --   stepper numbers during moving with constant speed
%   sn_d   --   stepper numbers of deceleration
%   pf_i   --   initial pulse frequency
%   pf_m   --   maximum pulse frequency
%   s_u    --   unit steps, default value is 1
%
% varargout:
%   f_list   --  frequencies list
%   dt_list  --  time periods list
%   
% copyright (c) wulx, <gurdy.woo@gmail.com>
% last modified by wulx, 2013/10/16

% default settings
if nargin < 6, s_u = 1; end
sn = sn_a + sn_c + sn_d;
f_list = zeros(1, sn); % pre-allocation for better performance
dt_list = zeros(1, sn);

% unit time period
dt = s_u / pf_i;

dt_list(1:end) = dt;

% #1 acceleration --------------------------------------------------------%
% divide frequencies uniformly in accordance with frequency range (pf_i, pf_m) and step numbers (sn_a)
f_list(1:sn_a) = linspace(pf_i, pf_m, sn_a);

% #2 moving with constant speed ------------------------------------------%
f_list(sn_a + (1:sn_c)) = pf_m;

% #1 deceleration --------------------------------------------------------%
f_list(sn_a + sn_c + (1:sn_d)) = linspace(pf_m, pf_i, sn_d);
