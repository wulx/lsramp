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

sn = sn_a + sn_d - 1;
f_list = zeros(1, sn); % pre-allocation for better performance
dt_list = zeros(1, sn); % ditto

% speed up -------------------------------------------------------%
% a_f -- frequency acceleration
[a_f, v_0, v_m] = ramp_eq(pf_i, pf_m, sn_a, s_u);

% motion formula:
%   v1^2 - v0^2 = 2*a*S
% rising edges
v_a = [v_0 arrayfun(@(n) sqrt(v_0^2 + 2*a_f*n*s_u), 1:sn_a)];

% reference line
%dt_a = (v_a(2:end) - v_a(1:end-1)) / a_f;
%t_list = arrayfun(@(n) sum( dt_a(1:n) ), 0:n_a);
%plot(t_list, v_list, 'k-')

f_list(1:sn_a) = (v_a(1:end-1) + v_a(2:end)) / 2;
dt_list(1:sn_a-1) = s_u ./ f_list(1:sn_a-1);

% at max speed ---------------------------------------------------%
%f_list(sn_a) = pf_m;
dt_list(sn_a) = (sn_c + 2)*s_u / pf_m;

% speed down -----------------------------------------------------%
[d_f, u_0, u_m] = ramp_eq(pf_i, pf_m, sn_d, s_u);

% similar to v_a, falling edges
u_d = [arrayfun(@(i) sqrt(u_m^2 - 2*d_f*(i-1)*s_u), 1:sn_d) u_0];

f_list(sn_a+1:end) = (u_d(2:end-1) + u_d(3:end)) / 2;
dt_list(sn_a+1:end) = s_u ./ f_list(sn_a+1:end);

if sum( strcmp(method, {'round', 'fix'}) ) == 1 % round or fix
    f_list = feval(method, f_list);
    
    s_list = s_u * ones(1, sn);
    s_list(sn_a) = (sn_c + 2)*s_u;
    dt_list = s_list ./ f_list;
end

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
