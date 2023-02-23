function [O,itout]=save_shorelinesSS(O,S,TIME,COAST,CLIFF,WAVE,TRANSP,STRUC,NOUR,GROYNE)
% function [O]=save_shorelines(O,it,dt,tc,nt,timenum1,x_mc,y_mc,dSds,x_mc1,y_mc1,PHIc,QS,HSo,PHIo,tper,HS,PHI,dPHI,HSbr,PHIbr,dPHIbr,hbr,STRUC.x_hard,STRUC.y_hard,x_nour,y_nour,GROYNE.x,GROYNE.y); 
%
% INPUT:
%      O                                                    % Output structure (where data is added to)
%      it,dt,tc,nt,timenum1,timenum0,itout,storageinterval  % time parameters
%      x_mc,y_mc,dSds,...                                   % x,y (after coastline update) and coastline change
%      x_mc1,y_mc1,PHIc,QS,...                              % x,y (before coastline update) and corresponding QS
%      HSo,PHIo,tper,PHIf,...                               % offshore waves & lower-shoreface orientation
%      HS,PHI,dPHI,...                                      % nearshore waves at depth-of-closure
%      HSbr,PHIbr,dPHIbr,hbr,...                            % waves at point of breaking
%      STRUC                                                % structures
%          .x_hard
%          .y_hard
%      x_nour,y_nour,...                                    % nourishment locations
%      GROYNE                                               % groyne locations
%          .x
%          .y
% OUTPUT:
%      O                                        % Output structure with corresponding fields to the input variables
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

    it=TIME.it;
    dt=TIME.dt;
    tc=TIME.tc;
    nt=TIME.nt;
    timenum1=TIME.tnow;
    timenum0=TIME.timenum0;
    itout=TIME.itout;
    storageinterval=S.storageinterval;        % time parameters
    x_mc=COAST.x_mc;
    y_mc=COAST.y_mc;
    dSds=COAST.dSds_mc;                     % x,y (after coastline update) and coastline change
    x_mc1=COAST.x_mc1;
    y_mc1=COAST.y_mc1;
    PHIc=COAST.PHIc_mc;
    QS=TRANSP.QS_mc;                          % x,y (before coastline update) and corresponding QS
    HSo=WAVE.HSo_mc;
    PHIo=WAVE.PHIo_mc;
    tper=WAVE.TP;
    PHIf=COAST.PHIf;                          % offshore waves & lower-shoreface orientation
    HS=WAVE.HStdp_mc;
    PHI=WAVE.PHItdp_mc;
    dPHI=WAVE.dPHItdp_mc;                     % nearshore waves at depth-of-closure
    HSbr=WAVE.HSbr_mc;
    PHIbr=WAVE.PHIbr_mc;
    dPHIbr=WAVE.dPHIbr_mc;
    hbr=WAVE.hbr_mc;                          % waves at point of breaking
    x_hard=STRUC.x_hard;
    y_hard=STRUC.y_hard;                      % structures
    x_nour=NOUR.x_nour;
    y_nour=NOUR.y_nour;                       % nourishment locations
    x_groyne=GROYNE.x;
    y_groyne=GROYNE.y;
    x_mc_cliff=CLIFF.x_mc;
    y_mc_cliff=CLIFF.y_mc;
    x_cliff=CLIFF.x;
    y_cliff=CLIFF.y;
    dndtcliff=CLIFF.dndt;
    retreatcliff=CLIFF.retreat;
    Kcgrid=CLIFF.Kc_grid;
    HSoff=CLIFF.HSo;
    TPeak=CLIFF.TP;
    if isfield(CLIFF,'TWL') && ~isempty(CLIFF.TWL)
        TWLevel=CLIFF.TWL;
    end
    Omm=CLIFF.Om;
    RetreatVolTotal=CLIFF.Retreat_Vol_total;
    RetreatVolToe=CLIFF.Retreat_Vol_toe;
    DriftVolShore=CLIFF.Drift_Vol_shore;
    % BeachW=COAST.BeachWidth;
    if isfield(S,'SurfZoneWidth') && ~isempty(S.SurfZoneWidth)
        COAST.ForeshoreSlope = abs(CLIFF.Ej + S.d)./ (COAST.BeachWidth + S.SurfZoneWidth);
        COAST.Beta = COAST.ForeshoreSlope;    
    else
        COAST.BeachSlope = CLIFF.Ej ./ COAST.BeachWidth;
        COAST.Beta = COAST.BeachSlope;
    end
    
    if (timenum1-timenum0)>=itout*storageinterval

        %% INITIALIZE OUTPUT FIELDS
        if it==0 
            fieldnmO = {'it','dt','tc','nt','timenum',...
                        'n','x','y','distx','distQS','dSds','dndt'...
                        'n1','x1','y1','distx1','distQS1','PHIc','QS'...
                        'HSo','PHIo','tper','PHIf',...
                        'HS','PHI','dPHI',...
                        'HSbr','PHIbr','dPHIbr','hbr',...
                        'x_hard','y_hard','n_hard',...
                        'x_nour','y_nour','n_nour',...
                        'x_groyne','y_groyne','x_mc_cliff',...
                        'y_mc_cliff','x_cliff','y_cliff',...
                        'retreat_cliff','Kc_grid','HSo','TP','TWL',...
                        'Om','Retreat_Vol_Total',...
                        'Retreat_Vol_Toe','Drift_Vol_Shore','n_mc_cliff',...
                        'n_mc1_cliff','distx_cliff','distx1_cliff',...
                        'distQS_cliff','distQS1_cliff','dndt_cliff',...
                        'BeachWidth','tan_B','SL','eta_setup',...
                        'R2_percent_Stockdon','cliffposition','coastposition'};   
            for ff=1:length(fieldnmO)
                O.(fieldnmO{ff})=[];
            end
        end  

        nans=find(isnan(x_mc));
        nans1=find(isnan(x_mc1));
        n_mc=length(nans)+1;
        n_mc1=length(nans1)+1;

        nans=find(isnan(x_mc_cliff));
        nans1=find(isnan(x_cliff));
        n_mccliff=length(nans)+1;
        n_mc1cliff=length(nans1)+1;

        %% compute distance x,y and transport locations: coast-line
        dx=x_mc(2:end)-x_mc(1:end-1);dx(isnan(dx))=0;
        dy=y_mc(2:end)-y_mc(1:end-1);dy(isnan(dy))=0;
        distx = [0,cumsum((dx.^2+dy.^2).^0.5)];
        distQS = (distx(1:end-1)+distx(2:end))/2;
        dx=x_mc1(2:end)-x_mc1(1:end-1);dx(isnan(dx))=0;
        dy=y_mc1(2:end)-y_mc1(1:end-1);dy(isnan(dy))=0;
        distx1 = [0,cumsum((dx.^2+dy.^2).^0.5)];
        distQS1 = (distx1(1:end-1)+distx1(2:end))/2;
        
        %% compute distance x,y and retreat locations: cliff-line
        dx_cliff=x_mc_cliff(2:end)-y_mc_cliff(1:end-1);dx_cliff(isnan(dx_cliff))=0;
        dy_cliff=y_mc_cliff(2:end)-y_mc_cliff(1:end-1);dy_cliff(isnan(dy_cliff))=0;
        distxcliff = [0,cumsum((dx_cliff.^2+dy_cliff.^2).^0.5)];
        distQScliff = (distxcliff(1:end-1)+distxcliff(2:end))/2;
        dx_cliff=x_cliff(2:end)-x_cliff(1:end-1);dx_cliff(isnan(dx_cliff))=0;
        dy_cliff=y_cliff(2:end)-y_cliff(1:end-1);dy_cliff(isnan(dy_cliff))=0;
        distx1cliff = [0,cumsum((dx_cliff.^2+dy_cliff.^2).^0.5)];
        distQS1cliff = (distx1cliff(1:end-1)+distx1cliff(2:end))/2;

        % store data in structure (which is automatically resized if needed)   
        [O.it]        = addVECTOR(O.it, it);
        [O.dt]        = addVECTOR(O.dt, dt);
        [O.tc]        = addVECTOR(O.dt, tc);
        [O.nt]        = addVECTOR(O.nt, nt);
        [O.timenum]   = addVECTOR(O.timenum, timenum1);

        % parameters after coastline change
        [O.n]         = addVECTOR(O.n, n_mc');
        [O.x]         = addVECTOR(O.x, x_mc');
        [O.y]         = addVECTOR(O.y, y_mc');
        [O.distx]     = addVECTOR(O.distx, distx');
        [O.distQS]    = addVECTOR(O.distQS, distQS');
        [O.dSds]      = addVECTOR(O.dSds, COAST.dSds');
        [O.dndt]      = addVECTOR(O.dndt, COAST.dndt');

        % parameters before coastline change
        [O.n1]        = addVECTOR(O.n1, n_mc1');
        [O.x1]        = addVECTOR(O.x1, x_mc1');
        [O.y1]        = addVECTOR(O.y1, y_mc1');
        [O.distx1]    = addVECTOR(O.distx1, distx1');
        [O.distQS1]   = addVECTOR(O.distQS1, distQS1');
        [O.PHIc]      = addVECTOR(O.PHIc, PHIc');
        [O.QS]        = addVECTOR(O.QS, QS');

        % parameters after cliff-line change
        [O.n_mc_cliff]      = addVECTOR(O.n_mc_cliff, n_mccliff');
        [O.x_cliff]         = addVECTOR(O.x_cliff, x_mc_cliff');
        [O.y_cliff]         = addVECTOR(O.y_cliff, y_mc_cliff');
        [O.distx_cliff]     = addVECTOR(O.distx_cliff, distxcliff');
        [O.distQS_cliff]    = addVECTOR(O.distQS_cliff, distQScliff');

        
        % if isfield(S,'x_mc') && ~isempty(S.x_mc) ...
        %         && isfield(COAST,'BeachWidth') && ~isempty(COAST.BeachWidth) 
        % if ~isfield(CLIFF,'Pos')
        %     CLIFF.Pos = CLIFF.x(1);
        %     COAST.Pos = COAST.x(1);
        % end
            [O.cliffposition]        = addVECTOR(O.cliffposition,CLIFF.x);%xc=CLIFF.x
            [O.coastposition]        = addVECTOR(O.coastposition,COAST.x);
        % elseif isfield(CLIFF,'x_mc') && ~isempty(CLIFF.x_mc) ...
        %         && isfield(COAST,'BeachWidth') && ~isempty(COAST.BeachWidth) 
        %     [O.cliffposition]        = addVECTOR(O.cliffposition,0+COAST.BeachWidth(end));   %0.00 wrt shoreline
        %                                                                                      % (end) for last transect
        % end

        % parameters before cliff-line change
        [O.n_mc1_cliff]       = addVECTOR(O.n_mc1_cliff, n_mc1cliff');
        [O.x_mc_cliff]        = addVECTOR(O.x_mc_cliff, x_mc_cliff');
        [O.y_mc_cliff]        = addVECTOR(O.y_mc_cliff, y_mc_cliff');
        [O.distx1_cliff]      = addVECTOR(O.distx1_cliff, distx1cliff');
        [O.distQS1_cliff]     = addVECTOR(O.distQS1_cliff, distQS1cliff');
%         [O.PHIc_cliff]      = addVECTOR(O.PHIc_cliff, PHIc_cliff');
%         [O.QS_cliff]        = addVECTOR(O.QS_cliff, QS_cliff');
        [O.retreat_cliff]     = addVECTOR(O.retreat_cliff, retreatcliff');
        [O.dndt_cliff]        = addVECTOR(O.dndt_cliff, dndtcliff');
        [O.Retreat_Vol_Total] = addVECTOR(O.Retreat_Vol_Total, RetreatVolTotal');
        [O.Retreat_Vol_Toe]   = addVECTOR(O.Retreat_Vol_Toe, RetreatVolToe');
        [O.Drift_Vol_Shore]   = addVECTOR(O.Drift_Vol_Shore, DriftVolShore');
        [O.Kc_grid]           = addVECTOR(O.Kc_grid, Kcgrid');
        [O.HSo]               = addVECTOR(O.HSo, HSoff');
        [O.TP]                = addVECTOR(O.TP, TPeak');
        if isfield(CLIFF,'TWL') && ~isempty(CLIFF.TWL)
            [O.TWL]   = addVECTOR(O.TWL, TWLevel');
        end
        [O.Om]                    = addVECTOR(O.Om, Omm');
        [O.BeachWidth]            = addVECTOR(O.BeachWidth, COAST.BeachWidth');
        [O.tan_B]                 = addVECTOR(O.tan_B, COAST.Beta);
        [O.SL]                    = addVECTOR(O.SL, CLIFF.SL);
        [O.R2_percent_Stockdon]   = addVECTOR(O.R2_percent_Stockdon, CLIFF.R2_percent_Stockdon);
        [O.eta_setup]             = addVECTOR(O.eta_setup, CLIFF.eta_setup);

        % offshore waves
        [O.HSo]       = addVECTOR(O.HSo, HSo);
        [O.PHIo]      = addVECTOR(O.PHIo, PHIo);
        [O.tper]      = addVECTOR(O.tper, tper);
        [O.PHIf]      = addVECTOR(O.PHIf, PHIf);

        % nearshore waves at depth-of-closure                      
        [O.HS]        = addVECTOR(O.HS, HS');
        [O.PHI]       = addVECTOR(O.PHI, PHI');
        [O.dPHI]      = addVECTOR(O.dPHI, dPHI');

        % waves at point of breaking
        [O.HSbr]      = addVECTOR(O.HSbr, HSbr);
        [O.PHIbr]     = addVECTOR(O.PHIbr, PHIbr);
        [O.dPHIbr]    = addVECTOR(O.dPHIbr, dPHIbr);
        [O.hbr]       = addVECTOR(O.hbr, hbr);

        % structures
        [O.x_hard]    = addVECTOR(O.x_hard, x_hard');
        [O.y_hard]    = addVECTOR(O.y_hard, y_hard');
        [O.n_hard]    = addVECTOR(O.n_hard, length(y_hard));

        % nourishments
        [O.x_nour]        = addVECTOR(O.x_nour, x_nour');
        [O.y_nour]        = addVECTOR(O.y_nour, y_nour');
        [O.n_nour]        = addVECTOR(O.n_nour, length(y_nour));

        % groyne locations
        [O.x_groyne]      = addVECTOR(O.x_groyne, x_groyne(:));
        [O.y_groyne]      = addVECTOR(O.y_groyne, y_groyne(:));
        
        itout=itout+1;
    end
end

% function to extend the size of a matrix
function [c]=addVECTOR(a,b)
    b=b(:);
    S1 = size(a);
    S2 = length(b);
    if S2>S1(1)
        a = [a;nan(S2-S1(1),S1(2))];
    elseif S1(1)>S2
        b = [b;nan(S1(1)-S2,1)];
    end
    c = [a,b];
end