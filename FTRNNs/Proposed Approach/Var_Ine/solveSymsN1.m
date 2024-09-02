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

Lamb = 8.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- Neural Network For The System -%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Projection Operator %
for i=1:length(x)
    Px(i,1) = ActFunc(Lb(i),Ub(i),x(i));
end
% Finite-Time Stabilizer (Lagrange Multiplier y)%
y = Stabilizer(Lb,Ub,x);
% Finite-time sliding mode operator (Lagrange Multiplier z) %
z = SMOperator(Lb,Ub,x,sig);

%- Proposed Projection Neural Networks -%
% Using Stabilizer to reach the Feasible space (For Optimization probems with convex functions)%
%- Dynamic for Neural Network -%
dx = Lamb*(-Fx(Px) + Px - x + A'*z + y);

%- Dynamic for sliding surface -%
dsig = Lamb*(-A*Fx(Px) + A*Px - A*x + A*A'*z + A*y);

%{
% ------------------------------------- %
%- Bi-layer Neural Network (RNN for VI) % (Ref. An extended projection neural network for constrained optimization
% ------------------------------------- %
% Projection Operator (Activation  Function) %
Ax = x - (W(x) - A'*sig);
for i=1:length(x)
    Px1(i,1) = ActFunc(Lb(i),Ub(i),Ax(i));
end
%- Neural Network -%
dx = Lamb*(Px1 - x);
dsig = Lamb*(-A*x + b);
%}

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

%- Finite Time Stabilizer Function -%
function Outy = Stabilizer(Lbx, Ubx, xin)
    % Auxiliar variables %
    Outy = zeros(length(xin),1); xl = norm(xin - Lbx);  xh = norm(xin - Ubx);
    % Set gains %
    k1 = 100; k2 = 100;
    
    % Select Case %
    for i = 1:length(xin)
        if xin(i) < Lbx(i)          % x < l
            Outy(i,1) = -k1*((xin(i) - Lbx(i))/xl) - k2*(xin(i) - Lbx(i));
        elseif  xin(i) > Ubx(i)     % x > h
            Outy(i,1) = -k1*((xin(i) - Ubx(i))/xh) - k2*(xin(i) - Ubx(i));
        else
            Outy(i,1) = 0;
        end
    end
end

%- Finite Time Sliding Modes Operator Function -%
function y = SMOperator(Lbx, Ubx, xin, sig)
    k1 = 3.5; k2 = 40.5;  
    % Auxiliar variables %
    y = zeros(length(sig),1); 
    % Set gains %
    % k1 = 10.4; k2 = 12.8;
    % Check Bounds and Select Case %
    if xin > Ubx           
        y = zeros(length(sig),1);
    elseif xin < Lbx
        y = zeros(length(sig),1);
    else
        % y = -k1*(sign(sig)) -k2*sig - k3*(sig.^3);
        y = -k1*(sig/norm(sig)) -k2*(sig);     % l < x < h
    end
        
end
