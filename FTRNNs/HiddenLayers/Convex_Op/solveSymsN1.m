function [dx,dsig,Px,y] = solveSymsN1(x,sig)
%%%%%%%%%%%%%%%%%%%%%
%- Getting States  -%
%%%%%%%%%%%%%%%%%%%%%
% Satates %
% x = Var(1:Nvar);
% sig = Var(Nvar+1:end);

% Nonconvex function %
% f(x) = -3*(x1^2) + 2*x1*x2 + 6*x1 -2*x2 - exp(x1) + exp(x2 + 2)
% f1(x) = exp(x1+1) - exp(x1) + 3
% Constraints %
% x1 -x2 = 1
% -2 <= x1,x2 <= 2
% ---------------------------------- %
% Optimal Values: x1*= -1, x2* = -2  %
% ---------------------------------- %

%- Lower Bound -%
Lb = [-2,-2]; 
%- Upper bound -%
Ub = [2,2]; 

% Matrix For equality Constraints %
A = [1, -1];
b = 1;

% Gradient For Objective Function %
% For variable x1 %
gx1 = exp(x(1)+1) - exp(x(1));
% For variable x2 %
gx2 = 0;
% Output Gradient %
gx = @(x)[gx1;gx2]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- Neural Network For The System -%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------------------------------- %
%- Bi-layer Neural Network (RNN for VI) % (Ref. An extended projection neural network for constrained optimization
% ------------------------------------- %
% Projection Operator (Activation  Function) %
Ax = x - (gx(x) - A'*sig);
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

%- Finite Time Stabilizer Function -%
function Outy = Stabilizer(Lbx, Ubx, xin)
    % Auxiliar variables %
    Outy = zeros(length(xin),1); xl = norm(xin - Lbx);  xh = norm(xin - Ubx);
    % Set gains %
    k1 = 100; k2 = 100;
    
    %- Function 3 -% 
    % k1 = 100; k2 = 120.5;
    % with Initial Values: 
    % x0 = [3,3,3]; z0 = [-7,-7,-7]; 
    
    %- Function 4 -% 
    % k1 = 0.5; k2 = 80.5;
    % with Initial Values: 
    % x0 = [-3,-3,-3,-3]; z0 = [-7,-7,-7]; 
    % x0 = [4,4,4,4]; z0 = [7,7,7,7]; 
    
    % k1 = 120.5; k2 = 80.5;
    % x0 = [3;5;2;6]; z0 = [7,7,7,7];
    
    % k1 = 35.6; k2 = 30.9;

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
    k1 = 5.5; k2 = 40.5;  
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
