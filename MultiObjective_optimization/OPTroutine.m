%%%%%%%%%%%%%%%%%%
%% OPTroutine.m %%
%%%%%%%%%%%%%%%%%%
% Implements the optimization routine for NNC algorithm.
%
% Xsqp [OUT]    : The decision vector calculated
% FUN  [OUT]    : Objective vector Value
% Flag [OUT]    : See fmincon help
% Options [OUT] : See fmincon help
% X     [IN]    : Decision Variable (initial guess)
% Dat   [IN]    : Parameters defined in NNCparam.m

%%%%%%%%%%%%%%%%%%%
%% Main Function %%
%%%%%%%%%%%%%%%%%%%
function [Xsqp, FLAG] = OPTroutine(x,Dat)   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Approach using the fmincon tool %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %-- See fmincon help for more information --%
    % options = optimset('Algorithm','active-set','Display','off');
    % [Xsqp, FUN, FLAG, Options] = fmincon(@(x) SQP(x,Dat),x,A,B,Aeq,Beq,lb,ub,@(x)SQPnolcon(x,Dat),options);
    % Xsqp: Optimal values founded for the respective fucnction
    % FUN: Respective function evaluated in Xsqp
    % FLAG: State of the algorimth
    % Options: Parameters used by the algorithm to minimize the function %
    %{
    lb = Dat.FieldD(:,1);       % lower bounds for x
    ub = Dat.FieldD(:,2);       % upper bounds for x
    Aeq = [];                   % Equality constraints of the form (Aeq)*(x)=Beq
    Beq = [];                   % 
    A = [];                     % Inequality constraints of the form (A)*(x) <= 0
    B = [];
    options = optimset('Algorithm','active-set','Display','off');
    [Xsqp,~,FLAG, ~] = fmincon(@(x) SQP(x,Dat),x,A,B,Aeq,Beq,lb,ub,@(x)SQPnolcon(x,Dat),options);
    jo = 0;                     % Only to set a breakpoint
    %}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Approach using the neural networks %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Xsqp = NNfun(x,Dat);
    % Only to test %
    FLAG = 1;   
    
    %-- Just in case we get NAN solution --%
    for nvar = 1:size(x,2)
        if isnan(Xsqp(1,nvar))
            Xsqp = x;                   
            break;
        end
    end  

end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Here we call the neural network function -%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = NNfun(x,Dat)
    nobj = Dat.nobj;
    Lb = Dat.FieldD(:,1);                                           % lower bounds for x
    Ub = Dat.FieldD(:,2);                                           % upper bounds for x
	gx = Dat.mopData.gx;                                            % Inequality constraints
  	grIne = Dat.mopData.grIne;                                      % Gradient of inequality constraints
    
    if Dat.LookingAnchors == 1
    % When looking for anchors points, the non-linear contraints required by NNC are not included %    
        grf = Dat.mopData.gr{nobj};                                 % Gradient of objective function (1 or 2)
        it = 1000;                                                  % Number of iterations %
    else
    % We are looking for Pareto points with the utopia plane constraint %    
        it = 4000;                                                  % Number of iterations %
        Auxgrf = Dat.mopData.gr{nobj};                              % Only auxiliar variable
        normgrf = @(x) Auxgrf(x)./Dat.MatrixL(nobj);                % Normalized gradient of objective function (1 or 2)
        grf = normgrf;                                             	% Gradient of objective function (1 or 2)
        
        %grf = Dat.mopData.gr{nobj};
        
        %- Additional constraints imposed by the NNC method -% 
        % normJ = (J-Dat.Utopia)./(Dat.MatrixL);
        % C = Dat.Nk*(normJ - Dat.normXpj)';                        % This one is required in NNC.
        
        J = Dat.mopData.func;                                       % Objective functions
        normJ = @(x)(J(x)-Dat.Utopia)./(Dat.MatrixL);               % Normalized objective functions 
        C = @(x) Dat.Nk*(normJ(x)- Dat.normXpj)';       % <= 0
        gx = @(x) [gx(x);C(x)];                                     % Extended vector of inequality constraints
        % gx = @(x) C(x);
                                             
        %-Gradient of the additional constraints imposed by the NNC method-%
        grf1 = Dat.mopData.gr{1};
        grf2 = Dat.mopData.gr{2};
        grC = @(x) (grf1(x)./(Dat.MatrixL(1)) - grf2(x)./(Dat.MatrixL(2)));
        grIne = @(x) [grIne(x), grC(x)];                            % Extended vector of gradients for inequality constraints
        % grIne = @(x) grC(x);  
    end
    % Time of integration %
    Ts = 50e-3;          
    % Slack variable for inequality constraints %
    y = rand(length(gx(x)),1);                                        
    % Slack variable for decision variables % 
    x = x';                                     
    for i=1:it
        %- Auxiliar variable (Note that we are not uding equality constraints in this case)-%
        xi = x - grf(x) - grIne(x)*y;   % A'*zs(i,:);         
        %- Projection Operator -%
        Px = zeros(length(xi),1);
        for n = 1:length(xi)
            Px(n,1) = piecewisex(Lb(n),Ub(n),xi(n));
        end
        % ----------------------------------------------------- %
        % Structure of Neural Networks (Differential Equations) %
        % ----------------------------------------------------- %
        %- Continuos- %
        dx = -x + Px;
        dy = -y + max((y + gx(x)),0);
        % PNNz(i,:) = -A*x(i,:)' + b;
        % ----------------------------------------------------- % 
        % Continouos Integration of the differential equations -%
        % ----------------------------------------------------- % 
        tspan = [0 Ts];
        [~,xout] = ode23(@(tx,xout) dx, tspan, x);
        [~,yout] = ode23(@(ty,yout) dy, tspan, y);
        x = xout(end,:)';
        y = yout(end,:)';
        % In this case we are not using equality constraints %
        % [tz,zout] = ode23(@(tz,zout) dz, tspan, z);
        % z = zout(end,:);
        
        % -------------------- %
        % Discrete-time Method %
        % -------------------- %
        %{
        dx = x + Ts*(-x + Px);
        dy = y + Ts*(-y + max((y + gx(x)),0 ));
        % dz = z + Tau*(-A*x + b);  
        %- Function Delay -%
        [x,y] = delayfnc(dx,dy);    
        %}
    end
    x = x';
end
%- Auxliar function to evaluate the projection -%
function y = piecewisex(Lbx, Ubx, xin)
    y = 0;
    idx1 = xin <= Lbx;
    y(idx1) = Lbx;
    idx2 = xin >= Ubx;
    y(idx2) = Ubx;
    y(~(idx1 | idx2)) = xin;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The non-linear constraints %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [C, Ceq] = SQPnolcon(x,Dat)
    J = Dat.mopData.func(x);
    % When looking for anchors points, the non-linear contraints required by NNC are not included %
    if Dat.LookingAnchors == 1      
        if Dat.NOBJ == 2
            C = [];
            Ceq = [];
        elseif Dat.NOBJ == 3
            C = [];
            Ceq = [];
        end
    % We are looking for Pareto points with the utopia plane constraint %        
    else                                                                                            
        if Dat.NOBJ == 2
            normJ = (J-Dat.Utopia)./(Dat.MatrixL);
            C = Dat.Nk*(normJ - Dat.normXpj)';                % This one is required in NNC.
            Ceq = [];
        elseif Dat.NOBJ == 3
            normJ = (J-Dat.Utopia)./(Dat.MatrixL);
            Dat.Nk;
            Dat.normXpj;
            C = Dat.Nk*(normJ-Dat.normXpj)';                %This one is required in NNC
            Ceq = [];
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Here we call the weighted function %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Jw = SQP(x,Dat)
    mop = Dat.mopData.func;
    J = mop(x);
    Jw = J*Dat.Mu';
end 