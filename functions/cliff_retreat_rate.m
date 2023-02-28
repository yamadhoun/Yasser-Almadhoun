function [CLIFF]=cliff_retreat_rate(S,COAST,CLIFF,WAVE,TIME)
% function [CLIFF]=cliff_retreat_rate(S,COAST,CLIFF,WAVE,TIME)
%
% TWL Kc=0.510; Hackney Kc=0.012; Trenhaile Kc=0.02/2.50; Energy_Flux Kc=0.20/0.10; Zhang Kc=0.035
%               
%% Copyright notice
%   --------------------------------------------------------------------
%   Copyright (C) 2020 IHE Delft & Deltares
%       
%       Dano Roelvink
%       d.roelvink@un-ihe.org
%       Westvest 7
%       2611AX Delft
%       
%       Bas Huisman
%       bas.huisman@deltares.nl
%       Boussinesqweg 1
%       2629HV Delft
%   
%   This library is free software: you can redistribute it and/or
%   modify it under the terms of the GNU Lesser General Public
%   License as published by the Free Software Foundation, either
%   version 2.1 of the License, or (at your option) any later version.
%   
%   This library is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   Lesser General Public License for more details.
%   
%   You should have received a copy of the GNU Lesser General Public
%   License along with this library. If not, see <http://www.gnu.org/licenses
%   --------------------------------------------------------------------
    
    %% INTERPOLATE OVER SPACE ADN TIME ALONG THE GRID
    
    % Get the interpolation method
    if ~isfield(S,'interpolationmethod') || isempty(S.interpolationmethod)
        disp('  Please set the interpolation method [S.interpolationmethod] in the main script ShorelineS...');
        return;
    end

    % Get the waves locations
    if isfield(S,'WVCfile') && ~isempty(S.WVCfile)
        WVC=WAVE.WVC; WC=WAVE.WC;
        xw=[]; yw=[];
        for kk=1:length(WVC)
            if (WAVE.spacevaryingwave)
               xw(kk,1)=WVC(kk).x;
               yw(kk,1)=WVC(kk).y; 
            end
        end
    end

    %% INTERPOLATE CALIBRATION COEFFICIENTS AND WATER DEPTH ALONG THE GRID
    
    % start of section

    % CALIBRATION COEFFICIENTS ALONG THE GRID: var1.Kc & var1.Kz
    % for test model and case study : Kc & d
    if isfield(S,'cliffModelCalibration') && ~isempty(S.cliffModelCalibration)
        if isfield(CLIFF,'Kc') && ~isempty(CLIFF.Kc)
            var1.Kc = CLIFF.Kc';
        else
            disp('  Please set the Kc calibration table [S.cliffModelCalibration] in the main script ShorelineS...');
            return;
        end
    elseif isfield(CLIFF,'Kconst') && ~isempty(CLIFF.Kconst) && ~isfield(S,'cliffModelCalibration')
        var1.Kc = CLIFF.Kconst *ones(length(COAST.xq));
    elseif isempty(CLIFF.Kconst) && ~isfield(S,'cliffModelCalibration')
        disp('  Please set the calibration constant [Kconst] in the text file paramsx...');
        return;
    else
        disp('  Please set the Kc calibration table [S.cliffModelCalibration] in the main script ShorelineSs,')
        disp('  or the calibration constants [Kconst] in the text file paramsx...')
        return;
    end
    
    % WATER DEPTH ALONG THE GRID: var1.d
    if isfield(CLIFF,'WaterDepth') && ~isempty(CLIFF.WaterDepth)
        for i=1:length(CLIFF.WaterDepth)
            if CLIFF.WaterDepth(i)<0    % or; CLIFF.WaterDepth(i)==-999 "-ve values"
                CLIFF.WaterDepth(i) = nan;
            end
        end
        CLIFF.WaterDepth = fillmissing(CLIFF.WaterDepth,'linear');
        var1.d = CLIFF.WaterDepth';
    end
    
    % var2. ALONG THE GRID: var2.useless
    if isfield(var1,'Kc') && ~isempty(var1.Kc)
        var2.useless = zeros(size(var1.Kc));
    end
    
    % CALIBRATION COEFFICIENTS AND WATER DEPTH ALONG THE GRID: Kc & d
    if isfield(S,'WVCfile') && ~isempty(S.WVCfile)
        [var1i,~,~] = get_interpolation_on_grid(S.interpolationmethod,CLIFF.x_mc,CLIFF.y_mc,xw,yw,var1,var2);
        if isfield(var1i,'Kc') && ~isempty(var1i.Kc)
            CLIFF.Kc_grid = var1i.Kc .* ones(size(CLIFF.x_mc));    
        end
        if isfield(var1i,'d') && ~isempty(var1i.d)
            CLIFF.d_grid = var1i.d;
        end
    elseif isfield(S,'Hso') && ~isempty(S.Hso) && isfield(S,'tper') && ~isempty(S.tper)
        CLIFF.Kc_grid=CLIFF.Kconst .* ones(size(CLIFF.x_mc));
    end

    %% INTERPOLAT WAVES IN SPACE FROM COASTLINE NEARSHORE TO CLIFF WALL POSITION
    %  INTERPOLATE THE WAVES CONDITIONS ALONG THE GRID: HSo & TP
    var3.HSo=WAVE.HSo;
    var3.TP=WAVE.TP;
    var4.useless=zeros(size(WAVE.TP));
    [var3i,~,~] = get_interpolation_on_grid(S.interpolationmethod,CLIFF.x_mc,CLIFF.y_mc,COAST.xq,COAST.yq,var3,var4);
    CLIFF.HSo=[]; CLIFF.TP=[];
    CLIFF.HSo=var3i.HSo; CLIFF.TP=var3i.TP;

    % end of section

    %% SWL = MSL + variational tide (waves) 
    if isfield(CLIFF,'MSL_times') && ~isempty(CLIFF.MSL_times)
        if isfield(CLIFF,'MSL') && ~isempty(CLIFF.MSL)
            SLevel = interp1(CLIFF.MSL_times, CLIFF.MSL, TIME.tnow);
        end
    end
    SL = SLevel * ones(size(CLIFF.x_mc));
    
    %% Beach slope: Beta
    if ~isfield(COAST,'BeachWidth') || isempty(COAST.BeachWidth) 
        COAST.BeachWidth = CLIFF.BaseBeach;
        CLIFF.Beta = abs(CLIFF.Ej) ./ (COAST.BeachWidth(1,1));
        CLIFF.Beta = CLIFF.Beta * ones(size(CLIFF.HSo));
    else
        CLIFF.Beta=[];
        CLIFF.Beta=COAST.BeachSlope';
    end

    %% Determine total water level: TWL
    Lw = 9.81/2/pi *CLIFF.TP.^2;
    R2_percent_Stockdon = 1.1 *0.35 .*CLIFF.Beta .*sqrt(CLIFF.HSo .*Lw) + ...
                          0.5 .*sqrt(CLIFF.HSo .*Lw .*(0.563 .*CLIFF.Beta.^2 + 0.004));
    eta_setup = 0.18 .*CLIFF.Beta .*CLIFF.TP .*sqrt(CLIFF.HSo *CLIFF.g);
    TWL = eta_setup + R2_percent_Stockdon + SL;                                  % SL = MSL (above NAVD88) + tidal variation (above MSL) + setup + R2%  + SLR

    %% Store water level data
    CLIFF.SL = SL;
    CLIFF.eta_setup = eta_setup;
    CLIFF.R2_percent_Stockdon = R2_percent_Stockdon;
    
    %% Nearshore depth of water
    if isfield(CLIFF,'dnearshore') && ~isempty(CLIFF.dnearshore)
        CLIFF.dnearshore(length(CLIFF.x_mc)+1:end) = [];
    elseif isfield(S,'dnearshore') && ~isempty(S.dnearshore)
        CLIFF.dnearshore = S.dnearshore * ones(1,length(CLIFF.x_mc));
    else
        CLIFF.dnearshore = 10 * ones(1,length(CLIFF.x_mc)); % typical nearshore depth=10m
    end

    %% Cliff erosion models: five formulae
    switch CLIFF.form
        case 1  % 'Zhang'
            C = Lw ./ CLIFF.TP;
            k = (2*pi)./Lw;
            Cg = C./2 .*(1+(2 .*k .*CLIFF.dnearshore) ./sinh(2.*k .*CLIFF.dnearshore));
            Fw = sind(CLIFF.theta)*(CLIFF.rho .*CLIFF.g .*CLIFF.HSo.^2 ./8) .*Cg;
            Fr = CLIFF.sigma .* CLIFF.Bc;
            for i = 1 : length(Fw)
                if Fw(i) >= Fr
                    dndt_cliff(1,i) = CLIFF.Kc_grid(i) *log(Fw(i)/Fr);
                else
                    dndt_cliff(1,i) = 0;
                end
            end
            CLIFF.TWL = TWL;
            CLIFF.Om = nan * ones(size(TWL));
            CLIFF.dndt = dndt_cliff;

        case 2  % 'TWL'
            CLIFF.TWL = TWL;
            CLIFF.Om = TWL / CLIFF.Ej;
            CLIFF.dndt = CLIFF.Om .* CLIFF.Kc_grid;

        case 3  % 'Trenhaile'        
            db = CLIFF.gamma .* CLIFF.HSo;
            Hb = Lw .*(0.142 *tanh((2*pi *db) ./(Lw)));                           % eqn name
            error = 100;                                                          % per cent or fraction
            nloops = 1;
            while max(error) > 20
                w = 1/CLIFF.alpha *(db + CLIFF.Ej + (CLIFF.retreat *CLIFF.alpha) - (TWL + SL));
                Om = CLIFF.rho *Hb / (2*CLIFF.gamma) .*exp(-CLIFF.Xdecay .*w) / 1000;
                dndt_cliff = Om .* CLIFF.Kc_grid;
                if nloops==1
                    dndt_cliff_prev = CLIFF.retreat;
                    error = abs(( dndt_cliff - dndt_cliff_prev ) ./ dndt_cliff_prev) * 100;
                else
                    dndt_cliff_prev = CLIFF.retreat - dndt_cliff;
                    error = abs(( dndt_cliff - dndt_cliff_prev ) ./ dndt_cliff_prev) * 100;
                end
                nloops = 1 + nloops;
            end
            CLIFF.TWL = TWL;
            CLIFF.Om = Om;
            CLIFF.dndt = dndt_cliff;

        case 4  % 'Energy_Flux'
            db = CLIFF.gamma .* CLIFF.HSo;
            Hb = Lw .*(0.142 *tanh((2*pi *db) ./(Lw)));                           % eqn name
            Cg = sqrt(CLIFF.g *db);
            error = 100;                                                          % per cent or fraction
            nloops = 1;
            while max(error) > 20
                w = 1/CLIFF.alpha *(db + CLIFF.Ej + (CLIFF.retreat *CLIFF.alpha) - (TWL + SL));
                Om = (1/8 *CLIFF.g *CLIFF.rho *Hb.^2 .*Cg) .*exp(-CLIFF.Xdecay .*w)/1000;
                dndt_cliff = Om .* CLIFF.Kc_grid;
                if nloops==1
                    dndt_cliff_prev = CLIFF.retreat;
                    error = abs(( dndt_cliff - dndt_cliff_prev ) ./ dndt_cliff_prev) * 100;
                else
                    dndt_cliff_prev = CLIFF.retreat - dndt_cliff;
                    error = abs(( dndt_cliff - dndt_cliff_prev ) ./ dndt_cliff_prev) * 100;
                end
                nloops = 1 + nloops;
            end
            CLIFF.TWL = TWL;
            CLIFF.Om = Om;
            CLIFF.dndt = dndt_cliff;

        case 5  % 'Hackney'
            CLIFF.Om = 1/8 *CLIFF.rho *CLIFF.g *(CLIFF.HSo + SL).^2 /1000;
            CLIFF.dndt = CLIFF.Om .* CLIFF.Kc_grid;
    end
    
    CLIFF.Retreat_Vol_total = CLIFF.dndt * CLIFF.Hc * 1 ;                           % Retreat total volume rate m3 / m / yr: volumetric rate
    CLIFF.Retreat_Vol_toe   = CLIFF.dndt * CLIFF.Hc * CLIFF.fc * 1;                 % Retreat volume rate m3 / m / yr remaining at cliff toe
    CLIFF.Drift_Vol_shore   = CLIFF.dndt * CLIFF.Hc * (1 - CLIFF.fc) * 1;           % Drift volume rate m3 / m / yr added to the alongshore sedment drift/transport
    CLIFF.Cfelev0           = CLIFF.Ej;                                             % Cliff toe elevation m
    
end
