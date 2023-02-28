function [CLIFF]=prepare_cliff(S,COAST)
% function [CLIFF]=prepare_cliff(S,COAST)
% 
% INPUT
%    S
%
% OUTPUT
%    CLIFF
%       .x_cliff
%       .y_cliff
%       .x_Cf     
%       .y_Cf
%       .Cfelev0
%
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
    
    fprintf('  Prepare cliffs \n');
    CLIFF=struct;
    CLIFF.ds0=S.ds0;
    CLIFF.smoothfac=S.smoothfac;
    CLIFF.XYoffset=[0,0];
    CLIFF.x_mc=[];
    CLIFF.y_mc=[];
    CLIFF.x0=[];
    CLIFF.y0=[];  
    CLIFF.Cfelev0=[];
    
    %% Read paramsx txt file: cliff and coast proporties
    fileID = fopen(S.cliffparams);
    C = textscan(fileID,'%s %s %f32');
    fclose(fileID);
    for i=1:size(C{1},1)
       CLIFF.(string(C{1,1}(i))) = C{1,3}(i);
    end

    %% Read cliff baseline: coordinates
    if isfield(S,'LDBcliff') && ~isempty(S.LDBcliff)
        xycoordcliff  = readtable(S.LDBcliff);
        CLIFF.x_mc = xycoordcliff.Var1';
        CLIFF.y_mc = xycoordcliff.Var2';
    elseif isfield(S,'x_cliff') && ~isempty(S.x_cliff) && isfield(S,'y_cliff') && ~isempty(S.y_cliff)
        CLIFF.x_mc=S.x_cliff;
        CLIFF.y_mc=S.y_cliff;
    elseif isfield(CLIFF,'BeachWidth') && ~isempty(CLIFF.BeachWidth)
        % see next for Offset Curve and Base Beach Width
    end

    %% Determine base beach required to calculate base Beta
    if isfield(S,'x_cliff') && ~isempty(S.x_cliff) && isfield(S,'y_cliff') && ~isempty(S.y_cliff)
        ax=CLIFF.x_mc; ay=CLIFF.y_mc;
        bx=COAST.x_mc; by=COAST.y_mc;
        [d_min, ~] = p_poly_dist(ax,ay,bx,by);
        CLIFF.BaseBeach = d_min;
    elseif isfield(S,'LDBcliff') && ~isempty(S.LDBcliff)
        ax=CLIFF.x_mc; ay=CLIFF.y_mc;
        bx=COAST.x_mc; by=COAST.y_mc;
        [d_min, ~] = p_poly_dist(ax,ay,bx,by);
        CLIFF.BaseBeach = d_min;
    elseif isfield(CLIFF,'BeachWidth') && ~isempty(CLIFF.BeachWidth) && ~isfield(S,'x_cliff') && ~isfield(S,'LDBcliff')
        CLIFF.BaseBeach = CLIFF.BeachWidth;
        % make offset curve to get CLIFF.x_mc & CLIFF.y_mc
        CLIFF.x_mc=[]; CLIFF.y_mc=[];
        [CLIFF.x_mc, CLIFF.y_mc] = offsetCurve(COAST.x_mc, COAST.y_mc, -CLIFF.BaseBeach);
    end

    %% Read water depth nearshore versus stations
    if isfield(S,'cliffNearshoreDepth') && ~isempty(S.cliffNearshoreDepth)
        stncliff  = readtable(S.cliffNearshoreDepth);                      % water depth vs stn nearhore for cliffs
        CLIFF.stn = stncliff.Var1';
        CLIFF.dnearshore = stncliff.Var2';
    end
    
    %% Read MSL txt file: time series
    if isfield(S,'cliffMSLevel') && ~isempty(S.cliffMSLevel)
        SLNAVD88 = readtable(S.cliffMSLevel);                              % MSL + Tide (wrt NAVD88 datum)
        SWLtime = datenum(SLNAVD88.Var1)';
        SWLevel = SLNAVD88.Var2';
        SLRise.SWLength = length(SWLevel);
    end
    
    %% Calculate SLR quadratic curve: SLRquad function
    SLRise.timestart =datenum(S.reftime);% datenum(SWLtime(1));            % wave.t(1,1); or SWLtime(1); or S.reftime                 % start time of the SLR curvce in [eg year]
    SLRise.timeend = datenum(S.endofsimulation);%(SWLtime(end));           % wave.t(end,1); or SWLtime(end); S.endofsimulation        % end time of the SLR curvce in [eg year]
    SLRise.SLRstart = 0.00;                                                % SLR at the start time of the SLR curvce in [m]
    SLRise.SLRend = 0.05; % 0.05 by default                                % SLR at the end time of the SLR curvce in [m]
    [SLRquad] = SLR_quad(SLRise);   
    
    %% Determine total MSL: MSL + SLR
    CLIFF.MSL = SWLevel' + SLRquad(:,2);
    CLIFF.MSL = CLIFF.MSL';
    CLIFF.MSL_times = SWLtime;

    %% Read calibration txt file: Kc factors table
    % [Kc] = Kc_calibration(stn)
    if isfield(S,'cliffModelCalibration') && ~isempty(S.cliffModelCalibration)
        Kc_calibrate = readtable(S.cliffModelCalibration);
        % assign Kc calibration values to CLIFF.stn
        for i = 1 : length(CLIFF.stn)
            for j = 1 : length(Kc_calibrate.Hindcast_stn)
                if CLIFF.stn(i) == Kc_calibrate.Hindcast_stn(j)
                    switch CLIFF.form
                        case 1  %'Zhang'
                            CLIFF.Kc(1,i) = Kc_calibrate.Kc_Zhang(j);
                        case 2  %'TWL'
                            CLIFF.Kc(1,i) = Kc_calibrate.Kc_TWL(j);
                        case 3  %'Trenhaile'
                            CLIFF.Kc(1,i) = Kc_calibrate.Kc_Tren(j);
                        case 4  %'Energy_Flux'
                            CLIFF.Kc(1,i) = Kc_calibrate.Kc_EF(j);
                        case 5  %'Hackney'
                            CLIFF.Kc(1,i) = Kc_calibrate.Kc_Hack(j);
                    end
                end
            end
        end
    else
        if isfield(CLIFF,'Kconst') && ~isempty(CLIFF.Kconst)
            CLIFF.Kc = CLIFF.Kconst;
        end
    end

end