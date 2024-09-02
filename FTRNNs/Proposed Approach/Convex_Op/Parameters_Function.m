clear all; clc;

%- Solver configuration -%
Ts = 5e-4;

% Function nonconvex problem %
% f(x) = -3*(x1^2) + 2*x1*x2 + 6*x1 -2*x2 - exp(x1) + exp(x2 + 2)
% Constraints %
% x1 -x2 = 1
% -2 <= x1,x2 <= 2
% ---------------------------------- %
% Optimal Values: x1*= -1, x2* = -2  %
% ---------------------------------- %
Nvar = 2;   % Number of Descition Variables
ECons = 1;  % Number of equality Constraints 

