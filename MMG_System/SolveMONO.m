function [OutM1,OutM2,fl] = SolveMONO(MG1,MG2,hr)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Parameters For Objective Functions %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %- Objective functions in the m-th MG within the MMG -%
    % **************************************** %
    % First Objective Function: Operating Cost %
    % **************************************** %
    % f1 =  sum(Costwt + Costpv + Costbb + CostDG + Costgrid)
    % for each interval time t with t = 1,...,24.
    
    % ****************************************** %
    % Second Objective Function (Emissions Cost) %
    % ****************************************** %
    % f2 =  sum(sum(xc*Edg) + xc*yg*Pgrid)
    % for each interval time t with t = 1,...,24.
            
    % ************************************************ %
    % Equality Operational Constraint (Energy Balance) %
    % ************************************************ %
    % Pload = varphi*Psh(m-j) = Pnc + varphi*Prec(m-j) + PUg
    
    % ***************************************************** %
    % Inequality Constraints For Projection Neural Networks %
    % ***************************************************** %
    % Omega = max(lambda(f1-f1*),(1-lambda)(f2-f2*)) <- This is the scalarization by Tchebyshev (nonsmooth)    
    % f1 - f1* - Omega <= 0 
    % f2 - f2* - Omega <= 0

    % This Physical Paramters Of Uncontrollable DERs (PVs and WTs) Are Given In Each MG  %
    % Minwt = 0;            % Kw/h   Min Power Generation For Wind Turbine         
    % Minpv = 0;            % Kw/h   Min Power Generation For Photovoltaic Panel 
    % Maxwt = 150;      	% Kw/h   Max Power Generation For Wind Turbine         
    % Maxpv = 200;       	% Kw/h   Max Power Generation For Photovoltaic Panel

    % [~,~,f1M,f2M] = Funct(w,Ivec);  % Call Auxiliar Function 
    % funGA = @(x)[f1M(x);f2M(x)];    % For Genetic Algorithms Technique %
    % options = gaoptimset('InitialPopulation',[0.5,0.4, 0.3]);   % Options For GA Method %
    % options = gaoptimset('PlotFcns',@gaplotpareto);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Input variables for optimization %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LoadMi = Load
    % PpvMi = PV Power Generated, PwtMGi = Wind Turbine Power Generated, SOCM = State Of Charge in BESSi 
    LoadM1 = MG1(1); PpvM1 = MG1(2); PwtM1 = MG1(3); SOCM1 = MG1(4);        % Variables For Microgrid 1 
    LoadM2 = MG2(1); PpvM2 = MG2(2); PwtM2 = MG2(3); SOCM2 = MG2(4);        % Variables For Microgrid 2 

    AuxSOC = [SOCM1;SOCM2];             % State of charge for the MMG system %

    % Note: Flag empty indicates that the code inside the conditional 
    % will run once during the simulation % 

    PowerPV = [PpvM1;PpvM2];            % PV power generated in each MG 
    PowerWT = [PwtM1;PwtM2];            % WT power generated in each MG 
    LMMG = [LoadM1;LoadM2];             % Load profiles in each MG 
    Pren = PowerPV + PowerWT;           % Available renewable power in each MG

    Pnet = Pren - LMMG;               	% Calculate power in each MG
    % If Pnet > 0, There Exists a Surplus Of Energy In MGi %
    % If Pnet < 0, There Exists a Deficit Of Energy In MGi %
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Main Loop For Multi-Objective Optimization %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % fprintf(['\n  ***************************** \n || Optimization In Hour: ' num2str(hr) ' || \n ***************************** \n ']);
    nMG = 2;                            % Number of MGs in MMG system
    nVar = 3;                         	% Number Of desicion variables 
    xOpNN = zeros(nMG,nVar);            % Auxiliar variables to save optimal values
    
    % fprintf('---------------------------- \n | Update NNs Inital States |\n ---------------------------- \n ');
    % Ener = zeros(nMG,1);
    % SBB = zeros(nMG,1);
    % ReqE = zeros(nMG,1);
    
    %- Obtain optimal solutions for each MG -%
    for M = 1:nMG                                                
        %- Information About Surplus/Deficit Of Energy In The MG -%
       	if Pnet(M,1) > 0
           Ener(M,1) = "Surplus";
        else
           Ener(M,1) = "Deficit";
        end
        %- Solve the scalarized multi-objective problem using the RNNs -%
        % Equality constraints % 
        Ic = abs(Pnet(M,1))/2;
        % Iincr = 20;
        Iincr = abs(Pnet(M,1))/20;
        
        if Pnet(M,1) > 0    
            AqM = -[1,1,1];
            % Initial conditions for descion variables %
            % xS1 = -[unifrnd(Ic-Iincr,Ic);unifrnd(0,1);unifrnd(Ic,Ic+Iincr)];
            xS1 = -[Ic-Iincr;0.5;Ic+Iincr];
        else
            AqM = -[1,1,1];
             %xS1 = [unifrnd(Ic-Iincr,Ic);unifrnd(0,1);unifrnd(Ic-Iincr,Ic)];
             xS1 = [Ic-Iincr;0.5;Ic+Iincr];
        end     
        bqM = Pnet(M);       
        
        % disp(xS1)  
        % Obtain optimal solutions for desicion variables %
        [xOpNN(M,:),sigma] = funcNN(xS1,Pnet(M),AqM,bqM,PowerPV(M),PowerWT(M),hr,AuxSOC(M));
        
        %disp("------ Sigma -------")
        %disp(sigma)
        %disp("--------------")
        
        %- Information About SOC Of BESSs -%
        if xOpNN(M,1) > 0
           SBB(M,1) = "Discharging";
        else
           SBB(M,1) = "Charging";
        end
        %- Information About Surplus/Deficit Of Energy After Optimization -%
        if xOpNN(M,3) > 0
           ReqE(M,1) = "Import";
        else
           ReqE(M,1) = "Export";
        end
    end
  
    Pgen = [Pren(1)+sum(xOpNN(1,:));            % Power Generated In The MMG System 
            Pren(2)+sum(xOpNN(2,:))];          
    EMMG = Pgen-LMMG;                           % Error between load and power generated 

    %- Print Values -%
    fprintf(' -------------------------------------- \n | Info about optimization in the first layer | \n -------------------------------------- \n');
    T = table(SBB,AuxSOC,Pnet,Pren,Ener,xOpNN, ReqE, Pgen,LMMG,EMMG,'VariableNames',{'|BB State|','|SOC|','|Pnet|','|Ren|','|Energy|','|Pbb, PDG, Pgrid| ','|Status|','|Pgen|','|Load|','|Error|'},'RowName',{'|| MG1 ||','|| MG2 ||'}); 
    disp(T)

    %%%%%%%%%%%%%%%%%%%
    %% Output Values %%
    %%%%%%%%%%%%%%%%%%%
    OutM1 = xOpNN(1,:);
    OutM2 = xOpNN(2,:);

    fl = hr +1;
end

% %%%%%%%%%%%%%%%%%%%% %
%- Auxiliar Functions -%
% %%%%%%%%%%%%%%%%%%%% %

% ************************** %
%- Neural Network Function -%
% ************************** %
function [xPos,Outsig] = funcNN(xS1,AuxP,A,b,Ppv,Pwt,hr,SOC)
 	%- Max and Min Power Generation  -%    
    Maxdg = (1e3)*150;                      % Kw/h   Max Power Generation For DG (Disel Generator)     
    Mindg = 1;                            	% Kw/h   Min Power Generation For Diesel Generator 
    Maxgrid = (1e3)*1500;                  	% Kw/h   Max Power Generation For The Utility Grid  
    peak = [(4:8),(16:21)];                	% (Hr) Period Range of Peak Tariff Category    
    if (SOC <= 30)
        Minbb = -(1e3)*30;                	% Kw/h   Min Power Generation For Battery Bank
        Maxbb = -100;
        Mingrid = 100;                     	% Kw/h   Min Power Generation For The Main Grid 
    else
        Minbb = -(1e3)*30;
        Mingrid = -(1e3)*1500;
        if ismember(hr,peak)
            Maxbb = (1e3)*(17);              	% Kw/h   Max Power Generation For BESS (Battery Energy Storage System)         
        else
            Maxbb = (1e3)*30;
        end

    end
    UpB = [Maxbb;Maxdg;Maxgrid];          	% Upper Bound
    LowB = [Minbb;Mindg;Mingrid];         	% Lower Bound  
    % Minbb = -(1e3)*200;                   % Kw/h   Min Power Generation For Battery Bank         
    % Mindg = 0;                         	% Kw/h   Min Power Generation For Diesel Generator     
    % Mingrid = -(1e3)*1500;             	% Kw/h   Min Power Generation For The Main Grid   
    
    %- Parameters of RNNs -%
    w = [0.5, 0.5];               	% Weight vector of dimension (pxm) 
                                  	% m: number of objective functions 
                                   	% p: number of desired solutions in the objective space
                                   	% We are giving the same importance to both objective functions 
  	
    [pNN,~] = size(w);            	% Number of neural networks used to characterize the Pareto front
    
    %- Iterations working without the break-%
    ite = 2500;                     % Max number of iterations to integrate differential equations
       
    nVar = 3;                      	% Number of desicion variables
    Ts = 5e-3;                      % Solver Time 

    Lamb = 1;                      	% Rate Of Convergence For The Neural Network 
    
    %- Call Scalarized Objective Function f(x), g(x) and Their Derivatives -%
    [ScalFunc,gScalFunc] = ObjFunct(w,Ppv,Pwt,hr,AuxP,SOC);          	
   
    %- Create periods of time to change gains values %
    Period1 = (14:17);
    Period2 = (18:19);
    
    %%%%%%%%%%%%%
    % Main Loop %
    %%%%%%%%%%%%%
    %- Calculate dynamics of recurrent neural networks -%
        for j=1:pNN
                Px = zeros(nVar,1);
                x = xS1;
                % Structure proposed by Xia in "An extended projection
                % neural network for constrained optimization" 
                % z = (1e-1)*ones(1,1);                                         % We use only z related to the equality constraints 

                %- Gains used in activation functions -%
                if ismember(hr,Period1)                     % Time slot from 14 hrs to 17 hrs
                    k2 = 150;                               % k2 = 150 (working value)
                    k1 = 62.5;                              % k1 = 62.5 (workin value)
                    sig = -(0.7)*ones(1,1);                 % sig = -(0.5)*ones(1,1);
                    
                elseif ismember(hr,Period2)                 % Time slot from 18 hrs to 19 hrs
                    k2 = 200;
                    k1 = 680;
                    sig = -(2.5)*ones(1,1); 
                else                                        % Other time slots  
                    k2 = 100;
                    k1 = 65.5;
                    sig = -(0.5)*ones(1,1); 
                end
                
                % Structure of our approach %     
                for it=1:ite                                                  	% Number of iterations for the integration                 
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %- Structure of recurrent neural networks -%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Structure of our approach %
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Finite-Time Stabilizer (Lagrange Multiplier y)%
                    y = Stabilizer(UpB,LowB,k2,x);
                    % Finite-time sliding mode operator (Lagrange Multiplier z) %
                    z = SMOperator(UpB,LowB,k1,x,sig);
                    
                	for n = 1:length(x)
                        Px(n,1) = piecewisex(LowB(n),UpB(n),x(n));          	% Apply the projection 
                    end
                    %- Dynamic for Neural Network -%
                    dx = Lamb*(-gScalFunc(Px) + Px - x + A'*z + y);
                    %- Dynamic for sliding surface -%
                    dsig = Lamb*(-A*gScalFunc(Px) + A*Px - A*x + A*A'*z + A*y);
                    
                    %- Integration of differential equations related to neural networks -%
                    tspan = [0 Ts];
                    [~,xout] = ode23(@(tx,xout) dx, tspan, x);
                    [~,sigout] = ode23(@(tz,sigout) dsig, tspan, sig);
                    x = xout(end,:)';
                    sig = sigout(end,:);

                end
        end

        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Output of optimal solutions %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %- Output using the methodlogy proposed by Xia in "An extended projection..." -% 
    % xPos = x;
    %- Ouput using our approach -% 
    %for n = 1:length(x)
    %    Px(n,1) = piecewisex(LowB(n),UpB(n),x(n));          	% Apply the projection 
    %end 
    xPos = Px;
    Outsig = sig;
end

% ********************************************************************** %
%- Finite Time Stabilizer Function (Used to reach the bound constraints)-%
% ********************************************************************** %
function Outy = Stabilizer(Ubx,Lbx,k,eta)
    % Auxiliar variables %
    Outy = zeros(length(eta),1); 
    xl = norm(eta - Lbx);  
    xh = norm(eta - Ubx);
    % Set gains %
      
    % Select Case %
    for i = 1:length(eta)
        if eta(i) < Lbx(i)          % x < l
            Outy(i,1) = -k*((eta(i) - Lbx(i))/xl);
        elseif  eta(i) > Ubx(i)     % x > h
            Outy(i,1) = -k*((eta(i) - Ubx(i))/xh);
        else
            Outy(i,1) = 0;
        end
    end
end

% *************************************************************** %
%- Finite time sliding mode operator:                             % 
% (Used to reach the sliding surface inside the bound constraints %
% *************************************************************** %
function Outz = SMOperator(Ubx,Lbx,k,eta,sig)
  	% Auxiliar variables %
    % Outz = zeros(length(sig),1); 
    % k1 = 5.5;
    
    % Check Bounds and Select Case %
    if eta > Ubx           
        Outz = zeros(length(sig),1);
    elseif eta < Lbx
        Outz = zeros(length(sig),1);
    else
        Outz = -k*(sig/norm(sig));     % l < x < h
    end
        
end

% ********************************** %
%- Evaluation of PieceWise Function -%
% ********************************** %
function y = piecewisex(Lbx, Ubx, xin)
    y = 0;
    idx1 = xin <= Lbx;
    y(idx1) = Lbx;
    idx2 = xin >= Ubx;
    y(idx2) = Ubx;
    y(~(idx1 | idx2)) = xin;
end

% ************************************************************************** %
%- Obtain a scalarized function and its derivative by using f1(x) and f2(x) -%
% ************************************************************************** %
function [ScalFunc,gScalFunc] = ObjFunct(w,Ppv,Pwt,hr,Pnet,SOC)
    %- Costs of grid energy -%
    peak = [(8:11),(16:21)];                % (Hr) Period Range of Peak Tariff Category 
    plain = (12:15);                        % (Hr) Period Range of Peak Tariff Category
    % valley = [(1:7),(22:24)];          	% (Hr) Period Range of Peak Tariff Category
    if ismember(hr,peak)
        % Peak tariff %
        Cgpur = 0.21;           % 0.17   	% (USD/Kwh) Cost of Purchasing Energy From The Grid 
        Cgsell = 0.12;          % 0.12     	% (USD/kWh) Cost of Selling Energy To The Grid
    elseif ismember(hr,plain)
       	% Plain tariff %
        Cgpur = 0.14;                       % (USD/Kwh) Cost of Purchasing Energy From The Grid 
        Cgsell = 0.061;         % 0.051             % (USD/kWh) Cost of Selling Energy To The Grid
    else
        % Valley tariff % 
        Cgpur = 0.062;                      % (USD/Kwh) Cost of Purchasing Energy From The Grid
        Cgsell = 0.039;                     % (USD/kWh) Cost of Selling Energy To The Grid
    end
    %- Emission values of CO2 for the main grid -%
    co2g = 428.55;                  	% g/Kwh CO2 emissions,  Other value: co2g = 388.55;  
    Eco2 = 0.00287;                    	% USD/Kg Cost by Penalization by CO2 Emissions
    bgg = co2g*Eco2;                  	% Cost by CO2 emissions by the main grid
    % so2g = 1.4;                      	% g/Kwh (SO2) sulfur dioxide emissions
    % noxg = 1.3;                      	% g/Kwh NOx emissions
    %- Coeffcients Costs Of Non-Dispatchable DERs -%
    Cwt = 0.0296;                      	% USD/Kwh   Cost of energy produced by the wind turbine
    Cpv = 0.0096;                     	% USD/Kwh   Cost of energy produced by the PV panel
    %- Coeffcients Costs of Battery Bank -%    
    %%%%%%%%%%%%%%%
    % Linear Form %
    %%%%%%%%%%%%%%%
    % Cbb = 0.082;                     	% USD/Kwh   Battery Bank Cost Generation
    
    %%%%%%%%%%%%%%%%%%
    % Quadratic From %
    %%%%%%%%%%%%%%%%%%
    % Values Obtained from Paper: 
    % Multiagent-Based Optimal Microgrid Control Using Fully Distributed Diffusion Strategy.
    aBBc = 5e-9;                       	% Original Value: aBBc = 11e-11; 
    bBBc = 3e-3;                    	% Original Value: bBBc = 8e-5;
    cBBc = 20;      
    om = 4;                           	% Original Value: om = 3;
    Pbbmax = 450e3;                   	% Maximum Charge or Discharge rate of the battery
    % Operating cost of energy produced by the battery bank %
    Cstbb = @(x) aBBc*(x(1) + om*Pbbmax*(1-(SOC/100)))^2 + bBBc*(x(1) + om*Pbbmax*(1-(SOC/100))) + cBBc;
    
    %- Coefficients Cost Of Diesel Generator -%   
    %%%%%%%%%%%%%%%
    % Linear Form %
    %%%%%%%%%%%%%%%
    % Paramaters obtained from paper: 
    % ("Efficient Energy Management of Renewable Resources in Microgrids").  
    % Cdg = 0.024;                      % USD/Kwh   Diesel Generator Cost Generation          
    % Cfuel = 1.2;                   	% USD/L Unit price of diesel
    % a = 0.0001;                             
    % b = 0.0438;                             
    % c = 0.3;
    
   % Emissions and Enviromental Vaule Of Diesel Generator %
    % co2 = 232.01;                 	% g/Kwh Emissions of carbono dioxide  (CO2) 
    % so2 = 0.464;                    	% g/Kwh Emissions of sulfur dioxide   (SO2)   
    % nox = 4.331;                    	% g/Kwh Emissions of nitrogen dioxide (NOx)
    % co = 2.32;                        % g/Kwh Emissions of carbon oxide     (CO)    
    % Eco = 0.125;                    	% USD/kg Penalization by emissions of CO2
    % Eso2 = 0.75;                    	% USD/kg Penalization by emissions of SO2  
    % Enox = 1;                       	% USD/kg Penalization by emissions of NOx 
    %%%%%%%%%%%%%%%%%%
    % Quadratic form %
    %%%%%%%%%%%%%%%%%%
    % Operating and environmental costs values obtained from the paper: 
    % Primal dual interior point algorithm for...
    % Operating costs %
    aDEc = 0.0547;
    bDEc = 1.7362;
    cDEc = 3.2456;
    % Environmenatl values %
    aDEe = 28.1444;
    bDEe = 1.7228;
    cDEe = 0.0017;
    % Operating cost of energy produced by the diesel generator %
    CstDE = @(x) aDEc*(x(2)^2) + bDEc*x(2) + cDEc;
    % Environmental cost of energy produced by the diesel generator %
    EnvCstDE = @(x) Eco2*(aDEe*x(2)^2 + bDEe*x(2) + cDEe);
    %- Check If Energy Will Be Purchased/Sold From/To The Grid -%1
    if Pnet < 0   
        Cg = Cgpur;
    else
        Cg = Cgsell;
    end
    % ************************** % 
    %- First Objective Function -%
    % ************************** %
    % Objective function given in a quadratic form %
    f1 = @(x) Cstbb(x) + CstDE(x) + Cg*x(3) + Cwt*Pwt + Cpv*Ppv;
    % *************************** %
    %- Second Objective Function -%
    % *************************** %
    % Objective function given in a quadratic form %
    f2 = @(x) EnvCstDE(x) + bgg*x(3);
    
    % *************************************************************** %
    %- Scalarized objective function (Using the weighted sum method) -%
    % *************************************************************** %
    ScalFunc = @(x) w(1)*f1(x) + w(2)*f2(x);
    
    %- Gradient of scalarized objective optimization function %
    gx1f1 = @(x) 2*aBBc*(x(1) + om*Pbbmax*(1-(SOC/100))) + bBBc;    % Gradient of first objective respect to x1
    gx2f1 = @(x) 2*aDEc*x(2) + bDEc;                                % Gradient of first objective respect to x2
    gx3f1 = Cg;                                                     % Gradient of first objective respect to x3
    grf1 = @(x) [gx1f1(x);gx2f1(x);gx3f1];
    
    gx1f2 = 0;                                                      % Gradient of second objective respect to x1
    gx2f2 = @(x) Eco2*(2*aDEe*x(2) + bDEe);                         % Gradient of second Objective Respect to x2
    gx3f2 = bgg;                                                    % Gradient of Second Objective Respect to x3
    
    grf2 = @(x) [gx1f2;gx2f2(x);gx3f2];
    gScalFunc = @(x) w(1)*grf1(x) + w(2)*grf2(x);                                     
end