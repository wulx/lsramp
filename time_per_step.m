function [f_list, dt_list] = time_per_step(sn_a, sn_c, sn_d, pf_i, pf_m, s_u, method)
%TIME_PER_STEP time per step algorithm
%
% varargin:
%   sn_a   --   number of stairsteps of speeding up
%   sn_c   --   number of stairsteps at constant speed (max speed)
%   sn_d   --   number of stairsteps of speeding down
%   pf_i   --   initial pulse frequency
%   pf_m   --   maximum pulse frequency
%   s_u    --   steps per stairstep, default value is 1
%   method --   approximation method
%
% varargout:
%   f_list   --  frequencies list
%   dt_list  --  time periods list
  
% copyright (c) wulx, <gurdy.woo@gmail.com>
% last modified by wulx, 2013/10/31

% check args
if any(~isfinite([sn_a, sn_c, sn_d]) | ([sn_a, sn_c, sn_d]<0))
    error('the number of steps should be finite non-negative.')
end

if any(~isfinite([pf_i, pf_m]) | ([pf_i, pf_m]<0))
    error('pulse frequencies should be finite non-negative.')
end

if ~isinteger(uint8(s_u)) || (s_u<1)
    error('steps per stairstep should be 8-bit positive integer (1 - 255).')
end

% default settings
if nargin < 7, method = 'ideal'; end
if nargin < 6, s_u = 1; end

sn = sn_a + sn_d + 1; % number of stairsteps

% number of stairsteps
s_list = s_u * ones(1, sn);
s_list(sn_a+1) = sn_c * s_u;

% frequencies list (pre-alloction)
f_list = nan(1, sn);

% speed up -------------------------------------------------------%
% a_f -- frequency acceleration
[a_f, v_0, ~] = ramp_eq(pf_i, pf_m, sn_a, s_u);

% motion formula:
%   v1^2 - v0^2 = 2*a*S
% rising edges
v_a = [v_0 arrayfun(@(n) sqrt(v_0^2 + 2*a_f*n*s_u), 1:sn_a)];

f_list(1:sn_a) = (v_a(1:end-1) + v_a(2:end)) / 2;

% at max speed ---------------------------------------------------%
f_list(sn_a+1) = pf_m;

% speed down -----------------------------------------------------%
[d_f, u_0, u_m] = ramp_eq(pf_i, pf_m, sn_d, s_u);

% similar to v_a, falling edges
u_d = [arrayfun(@(i) sqrt(u_m^2 - 2*d_f*(i-1)*s_u), 1:sn_d) u_0];

f_list(sn_a+2:end) = (u_d(1:end-1) + u_d(2:end)) / 2;
% ----------------------------------------------------------------%

% step time sequences
dt_list = s_list ./ f_list;

if any( strcmp(method, {'round', 'fix'}) ) % round or fix
    f_list = feval(method, f_list);
    
    dt_list = s_list ./ f_list;
end

end

function varargout = ramp_eq(f_i, f_m, n, s)
%RAMP_EQ linear speed ramp equation
% varargin:
%   f_i  --  initial frequency
%   f_m  --  maximum frequency
%   n    --  number of steps
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
