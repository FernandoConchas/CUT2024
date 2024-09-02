clear all; clc;

%- Solver configuration -%
Ts = 5e-4;

% *************** %
% Case Function 3 %
% *************** %
% variational Inequality Problem with Equality Constraints %
% 
%           [-x1 + 2x2 + x3]        [ 2]
%  U(x) =   [-x1 + x2 + 2x3]    d = [ 1]
%           [ x1 - x2 +  x1]        [-1]
% 0< = x1,x2,x3 <= 4
% Constraints %
% x1 - x2 + 2x3 = 1
% ------------------------------------------ %
% Optimal Values: x1*= 1, x2* = 0,   x3* = 0 %
% ------------------------------------------ %

Nvar = 3;   % Number of Descition Variables
ECons = 1;  % Number of equality Constraints 
