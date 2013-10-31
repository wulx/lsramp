function [f_list, dt_list] = steps_per_time(sn_a, sn_c, sn_d, pf_i, pf_m, s_u, method)
%STEPS_PER_TIME steps per time algorithm
%
% varargin:
%   sn_a   --   number of steps of speeding up
%   sn_c   --   number of steps at constant speed (max speed)
%   sn_d   --   number of steps of speeding down
%   pf_i   --   initial pulse frequency
%   pf_m   --   maximum pulse frequency
%   s_u    --   steps per stairstep, default value is 1
%   method --   approximation method
%
% varargout:
%   f_list   --  frequencies list
%   dt_list  --  time steps list

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

sn = sn_a + sn_d  + 1;
f_list = zeros(1, sn); % pre-allocation for better performance
s_list = zeros(1, sn); 

% unit time period
dt = s_u / pf_i;

switch method
    case 'ideal'
        dt_list = dt * ones(1, sn);
        dt_list(sn_a+1) = sn_c * dt;
        
        % #1 acceleration
        % divide frequencies uniformly in accordance with frequency range (pf_i, pf_m) and step numbers (sn_a)
        f_list(1:sn_a) = linspace(pf_i, pf_m, sn_a);
        
        % #2 moving with constant speed
        f_list(sn_a+1) = pf_m;
        
        % #3 deceleration
        % set f_list(sn_a) to pf_m again, duplicated but harmless!
        f_list(sn_a+2:end) = linspace(pf_m, pf_i, sn_d); 
    case {'round', 'fix'}
        % select round method
        %round_to = @round; % function handle
        %if strcmp(method, 'fix')
        %    round_to = @fix;
        %end
        round_to = @(x) feval(method, x); % more concise
        
        % #1 acceleration
        % divide frequencies uniformly and round them to nearest integer
        f_list(1:sn_a) = round_to(linspace(pf_i, pf_m, sn_a));
        s_list(1:sn_a) = round_to(dt * f_list(1:sn_a)); % step or pulse numbers should be integer
        
        % #2 moving with constant speed
        
        f_list(sn_a+1) = round_to(pf_m); % max frequency also should be rounded forcely
        s_list(sn_a+1) = round_to(sn_c*pf_m*dt);
        
        % #3 deceleration
        f_list(sn_a+2:end) = round_to(linspace(pf_m, pf_i, sn_d));
        s_list(sn_a+2:end) = round_to(dt * f_list(sn_a+2:end));
        
        dt_list = s_list ./ f_list;
    otherwise
        disp(['unknown approximation method: ' method])
        disp('available methods: ideal(default), round and fix')
end

end