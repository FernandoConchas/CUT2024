%%%%%%%%%%%%%%%%%%%%%%%%%
%% Overall Description %%
%%%%%%%%%%%%%%%%%%%%%%%%%
% This script generates the required parameters to run the NNC optimization algorithm 
% This code implements the NNC algorithm for 2 and 3 objectives as described in:
% A. Messac, A. Ismail-Yahaya and C.A. Mattson. The normalized normal 
% constraint method for generating the Pareto frontier structural and
% multidisciplinary optimization Volume 25, Number 2 (2003), 86-98.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Variables regarding the optimization problem %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

Objective = "Korn";
% Objective = "Other";

switch Objective
    case "Korn"
        funcVar = CostFunction1(Objective);
        NNCDat.CostProblem = 'Korn';                                % Cost function used
        
     	NNCDat.NOBJ = 2;                                            % Number of objectives
        NNCDat.NVAR = 2;                                            % Number of decision variables
        Lbound = [0;0];                                             % Lower bound
        Ubound = [5;3];                                             % Upper bound
        % Optimization bounds %
        NNCDat.FieldD = [ones(NNCDat.NVAR,1).*Lbound, ones(NNCDat.NVAR,1).*Ubound];              
        NNCDat.NRES = 0;                                            % Number of constraints
        % Important variables of objective function %
        NNCDat.mopData = funcVar; 
    otherwise
        disp("Something was wrong...")
        return
end
% Original Line of code %
% NNCDat.mop = str2func('CostFunction');            % Cost function file. In this
                                                    % file should be coded the
                                                    % calculation of the objective
                                                    % vector given a decision
                                                    % variable vector.
        
% Default values %
% NNCDat.CostProblem = 'DTLZ2';                     % Cost function used
% NNCDat.NVAR   = 10;                               % Number of decision variables
% NNCDat.FieldD =[zeros(NNCDat.NVAR,1), ones(NNCDat.NVAR,1)];     % Optimization bounds

% NNCDat.CostProblem = 'Korn';                        % Cost function used
% NNCDat.NVAR   = 2;                                  % Number of decision variables
% NNCDat.FieldD = [zeros(NNCDat.NVAR,1), [5;3]];       % Optimization bounds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Variables regarding the optimization algorithm %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NNCDat.mopOPT = str2func('OPTroutine');             % Optimization routine used. This
                                                    % call the optimization routine
                                                    % executed each time by the NNC
                                                    % algorithm. The fmincon routine
                                                    % is used in this application.

NNCDat.InitialGuess = [];                           % Use this if you have an initial
                                                    % guess to be used in the
                                                    % optimization. If empty, this
                                                    % implementation will use a
                                                    % random guess inside the
                                                    % optimization bounds.

%NNCDat.AnchorF=[1 0 1; 0 1 1; 1 1 0];              % Use this if you have a guess 
NNCDat.AnchorF = [];                                % on the anchor points. If empty,
NNCDat.AnchorX = [];                                % the NNC algorithm will 
                                                    % calculate them. 
                                                    % For the 3 objectives case the 
                                                    % anchor calculation is a 
                                                    % sensible step. Try the guess 
                                                    % (commented above) to see the 
                                                    % diferences on the outcome.
                                          
NNCDat.CounterGEN = 0;                              % Generation counter (Do not change).
NNCDat.CounterFES = 0;                              % Function Evaluation counter (Do not change).
                                          
NNCDat.Card = 20;                                   % Number of divisions in the 
                                                    % utopian line/plane. 
                                                    % For 3 objective problem a 
                                                    % suggested value is 10 (due to
                                                    % computational burden).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Variables regarding the Pareto Set filtering %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NNCDat.DominanceFiltering = 'yes';                  % This is a basic Dominance filtering;
                                                    % Write 'yes' if you
                                                    % want to filter the Pareto front
                                                    % and write 'no' if you do not
                                                    % want to do it.
                                          
NNCDat.SmartFiltering = 'yes';                      % This is the smart filter
NNCDat.RateFilter = [0.1 0.1];                      % Proposed by Matsson et Al. in:

% Mattson CA, Muller AA, Messac A. Smart Pareto Filter: Obtaining a minimal
% representation of multi-objective design space. Engineering Optimization 
% 2004 36(6):721-740.

%- Only will be executed if NNCDat.DominanceFiltering == 'yes'; -%
NNCDat.SaveResults = 'yes';                         % Write 'yes' if you want to 
                                                    % save your results after the
                                                    % optimization process;
                                                    % otherwise, write 'no';

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Call To NNC Function %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
OUT = NNC(NNCDat);                                  % Run the algorithm