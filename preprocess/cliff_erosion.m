function [dndt_cliff,Om,TWL] = cliff_erosion(Hs,Tp,SL,pars)
%
%
% global location location.Kz_zhang
% [Kc] = Kc_calibration(stn);
% pars.Kc = Kc;
% 
% global Kc
% pars.Kz = Kz ;
% pars.Kc = Kc * ones(size(Hs));
%
% CLIFF_EROSION
%
% Cliff erosion rate according to various formulations

%% INPUT:
%----------------------------Wave properties-------------------------------
%      Hs           : wave significant height in [meters]
%      Tp           : wave peak period in [seconds]
%      SL           : sea/water level in [meters]
%
%----------------------------Seawater properties---------------------------
%      g            : gravitational acceleration in [m/s2]
%      roh          : density of seawater in [kg/m3]
%
%      d            : depth of water in [m]
%      gamma        : wave breaking fraction [-]
%
%----------------------------Beach/Cliff properties------------------------
%      alpha        : beach slope at nearshore bet 3-10m in [degree]
%      Beta         : beach slope at foreshore in [degree]
%      theta        : angle bet waves and cliff vertical axis in [degree]
%
%      Ej           : cliff foot elevation in [meters]
%      sigma        : compressive strength of cliff components in [kg/m2]
%      Bc           : non-dimensional constant determined by the adherence of cliff components in [-]
%      dn_cliff_obs : beach slope at nearshore in [degree]
%
%----------------------------Model properties------------------------------
%      form         : relates to the type of model for cliff erosion in [-]
%      Kc           : calibration coefficient for the model in [-]
%      Kz           : Zhang constant with physical dimensions of a velocity [m/s]
%      Xdecay       : decay constant for wave energy dissipation in surf zone [0.03-0.07]


%% OUTPUT:
%      dndt_cliff   : retreat rate of cliffs in [m/yr]
%      Om           : cliff retreat forcing or accumulation of all wave energy over a given time in [kj/m2]
%      TWL          : total water level in [m]


%% Load paramsx txt file: input file
% paramsx = readtable('paramsx.txt');
% roh    = paramsx.Var3(1,1);
% g      = paramsx.Var3(2,1);
% % d      = paramsx.Var3(3,1);
% global d Beach Hc fc;
% gamma  = paramsx.Var3(4,1);
% alpha  = paramsx.Var3(5,1);
% Beta   = paramsx.Var3(6,1);
% theta  = paramsx.Var3(7,1);
% Beach  = paramsx.Var3(8,1);
% Ej     = paramsx.Var3(9,1);
% Hc     = paramsx.Var3(10,1);
% sigma  = paramsx.Var3(11,1);
% Bc     = paramsx.Var3(12,1);
% fc     = paramsx.Var3(13,1);
% % form   = paramsx.Var3(14,1);
% Kc     = paramsx.Var3(15,1);
% Kz     = paramsx.Var3(16,1);
% Xdecay = paramsx.Var3(17,1);


%% Calculations
% Wave length and total water level
Lw = 9.81/2/pi *Tp.^2;
R2_percent_Stockdon = 1.1 *0.35 .*pars.Beta .*sqrt(Hs .*Lw) + ...
                      0.5 .*sqrt(Hs .*Lw .*(0.563 *pars.Beta.^2 + 0.004));
eta_setup = 0.18 *pars.Beta .*Tp .*sqrt(Hs *pars.g);
TWL = eta_setup + R2_percent_Stockdon + SL; % SL = MSL (above NAVD88) + ...
%                                                  tide (above MSL) + setup + R + SLR
% Cliff erosion models
switch pars.form
    case 'Zhang'
        C = Lw ./ Tp;
        k = (2*pi)./Lw;
        Cg = C./2 .*(1+(2 .*k .*pars.d) ./sinh(2.*k .*pars.d));
        Fw = sind(pars.theta)*(pars.roh .*pars.g .*Hs.^2 ./8) .*Cg;
        Fr = pars.sigma .* pars.Bc;
        for i = 1 : length(Fw)
            if Fw(i) >= Fr
                dndt_cliff(i) = pars.Kz *log(Fw(i)/Fr);
            else
                dndt_cliff(i) = 0;
            end
        end
        dndt_cliff = dndt_cliff';
        TWL = TWL;
        Om(1:size(TWL)) = nan;
        Om = Om';

    case 'TWL'
        Om = (eta_setup + R2_percent_Stockdon + SL) / pars.Ej;
%         Om  = mean(Omi);
        dndt_cliff = Om .* pars.Kc;

    case 'Trenhaile'        
        db = pars.gamma .* Hs;
        Hb = Lw .*(0.142 *tanh((2*pi *db) ./(Lw)));    % eqn name
        dx = 0.1 *ones(size(Hs));                      % improve the first approximate
        error = 100;                                   % per cent or fraction
        nloops = 1;
        while error > 1
        w = 1/pars.alpha .*(db + pars.Ej + (dx .*pars.alpha) - (TWL + SL));
%         w = w(:,1);
        Om = pars.roh *Hb / (2*pars.gamma) .*exp(-pars.Xdecay .*w) / 1000;
        dndt_cliff = Om .* pars.Kc;
        error = abs(( dndt_cliff - dx ) / dx) * 100;
        dx = dx + dndt_cliff;
        nloops = 1 + nloops;
        end

    case 'Energy_Flux'
        db = pars.gamma .* Hs;
        Hb = Lw .*(0.142 *tanh((2*pi *db) ./(Lw)));    % eqn name
        Cg = sqrt(pars.g*db);
        dx = 0.1 *ones(size(Hs));                      % improve the first approximate
        error = 100;                                   % per cent or fraction
        nloops = 1;
        while error > 1
        w = 1/pars.alpha .*(db + pars.Ej + (dx .*pars.alpha) - (TWL + SL));
        Om = (1/8 *pars.g *pars.roh *Hb.^2 .*Cg) .*exp(-pars.Xdecay .*w)/1000;
        dndt_cliff = Om .* pars.Kc;
        error = abs(( dndt_cliff - dx ) / dx) * 100;
        dx = dx + dndt_cliff;
        nloops = 1 + nloops;
        end

    case 'Hackney'
        Om = 1/8 *pars.roh *pars.g *(Hs + SL).^2 /1000;
        dndt_cliff = Om .* pars.Kc;
end

end
