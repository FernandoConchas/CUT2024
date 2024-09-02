function [dx,dsig,Px,y] = solveSymsN1(x,sig)
%%%%%%%%%%%%%%%%%%%%%
%- Getting States  -%
%%%%%%%%%%%%%%%%%%%%%

% variational Inequality Problem with Equality Constraints %
% 
%           [-x1 + 2x2 + x3]        [ 2]
%  U(x) =   [-x1 + x2 + 2x3]    d = [ 1]
%           [ x1 - x2 +  x3]        [-1]
% 0< = x1,x2,x3 <= 4
% Constraints %
% x1 - x2 + 2x3 = 1
% ------------------------------------------- %
% Optimal Values: x1*= 1, x2* = 0,   x3* = 0 %
% ------------------------------------------- %
%- Lower Bound -%
Lb = [0,0,0]; 
%- Upper bound -%
Ub = [4,4,4]; 
% Matrix For equality Constraints %
A = [1,-1,-2];
b = 1;
d = [2;1;-1];

% In this case we use W as output function %
Fx = @(x) [-x(1) + 2*x(2) + x(3);
          -x(1) + x(2) + 2*x(3); 
           x(1) - x(2) + x(3)] + d;

% ------------------------------------- %
%- Bi-layer Neural Network (RNN for VI) % (Ref. An extended projection neural network for constrained optimization
% ------------------------------------- %
% Projection Operator (Activation  Function) %
Ax = x - (Fx(x) - A'*sig);
for i=1:length(x)
    Px(i,1) = ActFunc(Lb(i),Ub(i),Ax(i));
end
%- Neural Network -%
dx = Px - x;
dsig = -A*x + b;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluation of PieceWise Function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%- Activation Function -%
function y = ActFunc(Lbx, Ubx, xin)
    y = 0;
    idx1 = xin <= Lbx;
    y(idx1) = Lbx;
    idx2 = xin >= Ubx;
    y(idx2) = Ubx;
    y(~(idx1 | idx2)) = xin;
end