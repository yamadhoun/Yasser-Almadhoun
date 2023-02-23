% global S
% close all;
% clear all;clc;
% addpath(genpath('Shorelines_20221117\'))
addpath(genpath('..\functions\'))
addpath(genpath('..\test1\'))
figure (11);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODEL INPUT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S.cliffparams         = 'test1\paramsx.txt';
S.cliffNearshoreDepth = 'test1\Stations_Depth_1976_2005.txt';
S.cliffMSLevel        = 'test1\MSL_TIDE_NAVD88_dailyMax.txt';
S.reftime             = '2022-01-01';                                                    % Reference time <- leave empty to use t=0
S.endofsimulation     = '2052-12-31';                                                    % End time of simulation
S.outputdir           = 'Output_test1';  
S.interpolationmethod = 'weighted_distance';
S.Hso                 = 4;
S.tper                = 15;
S.x_mc                = [0 0];
S.y_mc                = [0 1000];
S.x_cliff             = [50 50];
S.y_cliff             = [0 1000];
S.xlimits             = [0 100];                                                         % x-limits of plot area
S.ylimits             = [0 1000];                                                        % y-limits of plot area
S.phiw0               = 330;                                                             % deep water wave angle [°N]
S.spread              = 90;                                                              % wave spreading [°] (wave_dir from range:  S.phiw0 +/- 0.5*S.spread)
S.d                   = 10;                                                              % Active profile height [m]
S.ddeep               = 25;                                                              % Offshore water depth for refraction
S.dnearshore          = 12;                                                              % Nearshore water depth for refraction
S.trform              = 'CERC';                                                          % Switch for transport formulation (e.g. S.trform='CERC', 'KAMP', 'MILH' or 'VR14')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Adaptaive Time Step: ATS: S.b>0 & 1>S.tc>0
% S.tc                  = 0.5;                                                           % factor to modify automatic timestep
% S.b                   =.5e6;                                                           % CERC : coeff in simple cerc formula
% S.d50                 = 6.0e-4; 
%% No Adaptaive Time Step: No ATS: S.b=0 & S.tc=0
S.tc                  = 0;
S.b                   = 0;
S.dt                  = 1/365;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S.growth              = 0;                                                               % calibration of nourishment growth rate
S.nourish             = 0;                                                               % switch (0/1) for nourishments
S.LDBnourish          = '';                                                              % LDB with nourishment locations ([Nx2] ASCII FILE WITHOUT HEADER) (i.e. polygon around relevant grid cells) <- leave empty to use interactive mode!
S.nourrate            = 100;                                                             % nourishment rate in polygon (m3/m/y)
S.nourstart           = 0;                                                               % start time of nourishment (y)
S.spit_width          = 50;                                                              % width of tip of spit
S.twopoints           = 0;                                                               % upwind treatment involving two points (1) or 1 point (0)
S.storageinterval     = 365;                                                             % output directory for plots and animation
S.plotinterval        = 1;
S.fignryear           = 1;
S.smoothfac           = 0.1;
S.phiw0               = S.phiw0 *pi/180;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RUN SHORELINES MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[S,O]=ShorelineS(S);
