function [f_list, dt_list] = time_per_step(sn_a, sn_c, sn_d, pf_i, pf_m, s_u, method)
%TIME_PER_STEP time per step algorithm
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
%   
% copyright (c) wulx, <gurdy.woo@gmail.com>
% last modified by wulx, 2013/10/17

% default settings
if nargin < 7, method = 'ideal'; end
if nargin < 6, s_u = 1; end


end

function varargout = ramp_eq(f_i, f_m, n, s)
%RAMP_EQ linear speed ramp equation
% varargin:
%   f_i  --  initial frequency
%   f_m  --  maximum frequency
%   n    --  stepper numbers
%   s    --  unit steps
% varargout:
%   a    --  acceleration or deceleration
%   v_0  --  starting frequency of reference line
%   v_m  --  maximum frequency of reference line
% copyright (c) wulx, <gurdy.woo@gmail.com>
% last modified by wulx, 2013/10/17

a = sym('a', 'positive');
syms v_0 v_m real
% solve multiple equations of linear speed ramp
sols = solve(v_m == f_m + 0.5*a*s/f_m, ...
    v_0 == f_i - 0.5*a*s/f_i, ...
    v_m^2 - v_0^2 == 2*a*n*s, ...
    a, v_0, v_m);

res = structfun(@(s) double(s), sols);
varargout = {res(1), res(2), res(3)};
end