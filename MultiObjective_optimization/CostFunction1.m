%% CostFunction.m 
% J  [OUT] : The objective Vector. J is a matrix with as many rows as
%            trial vectors in X and as many columns as objectives.
% x   [IN] : Decision Variable Vector. X is a matrix with as many rows as
%            trial vector and as many columns as decision variables.
% Dat [IN] : Parameters defined in NNCparam.m

%%%%%%%%%%%%%%%%%%%%
%% Main call Body %%
%%%%%%%%%%%%%%%%%%%%
function J = CostFunction1(Case)
    if strcmp(Case,"Korn")
        J = Korn();
    elseif strcmp(Case,"Korn")
        J = Korn();                                         % Here comes the call for a cost 
    end      
end

%%%%%%%%%%%%%%%%%%%%%%
%% Custom functions %%
%%%%%%%%%%%%%%%%%%%%%%
% Korn Function %
% Optimal values: xs1* = (0,0), xs2* = (5,3)
function J = Korn()
    %- Objective function 1 -%
    f1 = @(x) 4*(x(1)^2) + 4*(x(2)^2);
    %- Objective function 2 -%
    f2 = @(x) (x(1)-5)^2 + (x(2)-5)^2;
    % Gradient of f1 %
    grf1 = @(x) [8*x(1);8*x(2)];
    % Gradient of f2 %
    grf2 = @(x) [2*(x(1)-5);2*(x(2)-5)];
    % Establish the inequality constraints %
    g1 = @(x) (x(1)-5)^2 + (x(2))^2 - 25;   % <= 0
    g2 = @(x) -(x(1)-8)^2 -(x(2)+3)^2 + 7.7;   % <= 0
    gx = @(x) [g1(x);g2(x)];
    
    % For variable x1 %
    gr1x1 = @(x) 2*(x(1)-5);
    gr2x1 = @(x) -2*(x(1)-8);    

    % For variable x2 %
    gr1x2 = @(x) 2*x(2);    
    gr2x2 = @(x) -2*(x(2)+3);

    % Auxliar variable that contains the gradient of all decision variables %
    grIne = @(x) [gr1x1(x), gr2x1(x);gr1x2(x), gr2x2(x)];

    % Make a struct with the necesary variables %
    J.func = @(x) [f1(x), f2(x)];
    J.gr{1} = grf1;
    J.gr{2} = grf2;
    J.gx = gx;
    J.grIne = grIne;
end