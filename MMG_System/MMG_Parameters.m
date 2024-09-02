clear all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solver Configuration %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
TimeS = 24;                              % [Sec] Time For Simulation
Ts = 8e-4;                              % [Sec] General Sample Step
% Ts_C = 10*Ts;                         % [Sec] Sample Step 
% Ts_Control = 2*Ts;                    % [Sec] Sample Step
% Ts_Power = ;                          % [Sec] Sample Step For Controller
Time = 0:1:23;                          % Time Slots For a day 

%%%%%%%%%%%%%%%%%%%%%
%% Parameters Grid %%
%%%%%%%%%%%%%%%%%%%%%
Vgrid = 25e3;                           % [V]       Voltage Grid
Vs = 120e3;                             % [V]       Main Grid Voltage 
fr = 60;                                % [Hz]      Grid frecuency Voltage 
wn = 2*pi*fr;                           % [rad/s]   Grid frecuency Voltage 
Np = 250e6;                             % [w]       Nominal power 
VMic = 575;                             % [V]       Microgrid Voltage
NPMic = 150e6;                          % [w]       Nominal Power For Microgrid Transformers
Lgs = 600e-6;                           % [H]       Grid side converter Inductance Simulation: 400e-6 (H)  %Real values 0.3495 

% Rgs = 20e-6;   % Original value       % [Ohm]     Grid resistence  Simulation: 20e-6 (ohm)   %Real values 0.108
Rgs = 20e-6;

MicroVol = 575;                         % [V]       Microgrids Voltage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters Wind Turbibe (WT) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CWTdc = 1100e-6;                      % Capacitor DC Link Wind Turbine

MaxWind = 18;                           % [m/s] Max Wind Velocity     
MinWind = 3;                            % [m/s] Min Wind Velocity
NomWind = 10;                           % [m/s] Nominal Wind Velocity

% PnomWT = 100e3;  % (Original value)   % [W]   Nominal Power Wind Turbine
PnomWT = 400e3;                         % [W]   Nominal Power Wind Turbine


% Wind Profile For Each MG Within The MMG System (Using Database MERRA) %
WindM1 = [2.23,2.3,2.47,2.67,2.93,3.08,3.19,3.64,3.44,3.36,2.32,1.58,1.34...
    1.23,1.04,0.95,0.87,1.2,0.65,0.46,1.1,0.78,1.09,1.4]+3.2; %+3.8

WindM2 = [2.73,2.75,2.75,2.76,2.88,3.15,3.88,3.68,4.34,4.57,4.46,4.28,4.11...
    3.94,3.75,3.45,2.94,2.28,2.16,2.21,2.36,2.51,2.66,2.77]+1.5; %+2.5

WindM3 = [2.87,2.56,2.23,1.91,1.69,1.57,1.99,1.56,1.43,1.92,2.34,2.37,2.24...
    2.08,2.01,1.9,1.86,1.59,1.54,1.56,1.33,1.81,2.1,2.29]+2; % +2.5

WindM4 = [2.78,2.75,2.69,2.79,3.07,2.92,2.43,2.86,1.94,1.05,1.26,1.82,1.77...
    1.48,1.25,1.17,1.05,0.78,0.7,0.95,1.43,1.88,2.17,2.3]+4.1; % +4.4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parmeters Battery Bank %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NomVol = 600;                           % [V]       Nominal Voltaje
Pnombb = 400;                           % [Kw]      Nominal Power

% Original values %
% RtCapM1 = (1.1)*1e-1;                   % [Kw/h]    Rated Capacity Microgrid 1
% RtCapM2 = (1.1)*1e-1;                   % [Kw/h]    Rated Capacity Microgrid 2
% RtCapM3 = (1.1)*1e-1;                   % [Kw/h]    Rated Capacity Microgrid 3
% RtCapM4 = (1.1)*1e-1;                   % [Kw/h]    Rated Capacity Microgrid 4

RtCapM1 = (5.1)*1e-1;                   % [Kw/h]    Rated Capacity Microgrid 1
RtCapM2 = (5.1)*1e-1;                   % [Kw/h]    Rated Capacity Microgrid 2

Ah = 80;                                % [A/h]     Battery Capacity 
SOC = [70,80,75,65];                    % [%]       Initial Battery SOC
BatRes = 30;                            % [Sec]     Battery Response Time

%SOC = [70,60,50,80];                    % [%]       Initial Battery SOC
% RtCapM1 = (9.5)*1e-2;   % (7.6)*1e-2;   % [Kw/h]    Rated Capacity Microgrid 1
% RtCapM2 = (1.3)*1e-1;   % (1.1)*1e-1;   % [Kw/h]    Rated Capacity Microgrid 2
% RtCapM3 = (1.7)*1e-1;   % (1.1)*1e-1;   % [Kw/h]    Rated Capacity Microgrid 3
% RtCapM4 = (9.0)*1e-2;   % (1.1)*1e-1;   % [Kw/h]    Rated Capacity Microgrid 4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parmeters Photovoltaic Panel %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data From Merra2 NASA Dataset % 
%- Average Irradiation (W/m^2) and Temperature (C) Per day in MG1 -%
IrMG1 = [0,0,0,0,0,0,8.4237,87.6102,230.5238,455.5221,593.5651,678.5624,705.2509...
    663.6362,477.459,271.3045,150.3972,72.6094,20.9706,0,0,0,0,0]; %+10

TempMG1 = [24.76,24.45,24.08,23.73,23.47,23.22,23.08,23.64,24.81,26.92,28.17,28.95...
    29.53,29.87,29.73,29.27,28.95,28.48,27.79,27.18,26.76,26.23,25.69,25.06];
    
%- Average Irradiation (W/m^2) and Temperature (C) Per day in MG2 -%    
IrMG2 = [0,0,0,0,0,13.1616,158.4647,384.4585,594.6108,766.0255,876.8519,905.9282,857.1017...
    748.1747,596.2939,405.4984,198.3632,41.9603,0,0,0,0,0,0]; % -30

TempMG2 = [24.67,24.34,24.05,23.81,23.59,23.47,24.45,26.33,28.46,29.79,30.63,31.12,31.35...
    31.33,31.15,30.78,29.99,28.43,27.02,26.5,25.89,25.31,24.84,24.47];

%- Average Irradiation (W/m^2) and Temperature (C) Per day in MG3 -% 
IrMG3 = [0,0,0,0,0,6.1254,85.4846,231.4892,354.5752,503.6502,627.6343,679.0825,607.7515...
    436.9741,320.5312,210.0829,90.2474,15.0025,0,0,0,0,0,0]-40;

TempMG3 = [25.78,25.55,25.35,25.16,25.03,24.99,25.86,27.49,28.78,29.58,30.3,30.93,31.18...
    30.82,30.34,29.85,29.14,28.03,27.24,26.86,26.66,26.74,26.54,25.95];

%- Average Irradiation (W/m^2) and Temperature (C) Per day in MG4 -% 
IrMG4 = [0,0,0,0,0,7.1224,124.4849,344.2311,551.9198,724.0964,840.6983,840.6735,820.9622...
    731.0239,580.2064,367.0017,160.7112,28.1318,0,0,0,0,0,0]-70;

TempMG4 = [23.77,23.52,23.35,23.3,23.22,23.17,24.37,26.35,28.96,30.34,30.99,31.32,31.58...
    31.63,31.47,31.11,30.87,29.79,28.7,27.69,26.74,25.99,25.41,24.95];

%%%%%%%%%%%%%%%%%%%%%%
%% Diesel Generator %%
%%%%%%%%%%%%%%%%%%%%%%
PnomDG = 150e3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loads profiles in The MMG %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- Loads According Dataset (NREL : oedi-data-lake/nrel-pds-building-stock/end-use-load-profiles-for-us-building-stock/2023/comstock_amy2018_release_2/timeseries_aggregates/by_state/upgrade%3D12/state%3DAL) 
% Small Office (kW/h) %
PLoad1 = (0.18*1e1)*[19351.0770,18581.4994,27208.8579,32105.0767,30321.8573,41278.4873,45229.3026,48860.6062,49899.5490,54610.0882,53063.2429...
    53237.0729,49718.2848,44841.6642,43924.4254,41204.1987,39532.7488,37347.7957,41576.9259,41107.9637,40618.5496,39909.9045,35109.9421,29952.7555];

% Small Office (kW/h) %
%PLoad2 = [13529.6285,13145.3203,13707.0246,14370.7285,14443.6721,14394.7880,13930.8194,13254.3859,18137.7475,18057.6554,18153.0698,18466.0310...
%    18753.4387,18917.8711,18877.5062,18775.6975,18546.7200,18481.7277,18619.5634,18713.7441,18945.0101,18824.1704,18723.8391,13923.8398];

% Medium Office (kW/h) %
PLoad2 = (0.3*1e1)*([30234.381,34589.24,41714.573,45958.26,43639.216,49196.168,53175.33,62512.369,60959.416,52304.677,68173.109,54096.721...
    45892.723,34835.464,21588.854,18044.279,16418.323,15063.06,16581.566,17643.543,18267.403,20050.844,25847.36,30962.048]);  % -5000 

% Medium Office CA(kW/h) %
%PLoad3 = [10589.31, 12540.83, 14305.72, 16262.38, 18240.57, 19869.95, 20796.76, 23796.96, 24649.95, 21371.66, 23981.63...
%    25949.53,21861.54,17835.75,13589.02,9009.08,10120.44,12016.05,9711.82,8509.68,9815.41,9218.01,10188.58,11313.8];

% Retail Strip Mall (kW/h)%
%PLoad4 = [58508.7320,61912.6841,65573.8071,68642.2285,71026.0334,72148.2917,75971.8673,72581.2924,69361.1112,63473.4277,65491.0429,58832.5098...
%    46285.5460,45873.3562,47090.8346,38464.9778,41096.4581,45437.7744,48168.2419,57081.3370,55075.2720,65137.8578,68660.7543,61205.6595]-18000;

% Small Hotel (kW/h)
%PLoad3 = [3121.7353,2719.6772,2346.4347,2059.5081,1925.9906,1900.5063,2513.9303,2833.4486,3232.4102,3512.0676,3714.1053,3927.0109,4019.6531...
%    4071.2381,4090.5607,4156.0713,4201.9013,4311.3154,4391.7063,4452.6817,4051.3188,3771.5288,3434.4852,3209.7506];


%%%%%%%%%%%
%% Clock %%
%%%%%%%%%%%
ClockS = 0:0.2:TimeS;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check Time Execution For The Model %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% open('C:\Users\ferna\OneDrive\Escritorio\Neural-Simulations\Working-IN\Multi-Micro\MMG_Working_4MGs\4MMG_2Test\simpleMicrogrid_R2019B.slx');
% disp('Model is still running!')
% tic;
% sim 'simpleMicrogrid_R2019B.slx'
% executionTime=toc
% disp('Model is finished.')

