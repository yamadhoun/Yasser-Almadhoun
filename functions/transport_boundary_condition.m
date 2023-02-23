function [TRANSP]=transport_boundary_condition(TRANSP,COAST,TIME,GROYNE)
% this function applied for non cyclic shorelines 
% In case of more than one (non-cyclic) section with different...
% ... boundary the code should be adjusted 
%
%'neumann'...Neumann boundary
%'CTAN'...Constant orientation 
%'func'...Dirichlet boundary 
%'closed'...Dirichlet (wall boundary)
%'periodic'...Periodic boundary Q+position   % testing phase 
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

    eps=1e-6;
    idgroyne=[];
    %% start point (left)
    if ~isempty(GROYNE.x)
       [idgroyne]=find(GROYNE.idcoast(:,2)==COAST.i_mc);
    end
    if isempty(idgroyne)  % Do not apply any of these transport boundary conditions 
                        % in case of a groyne on this side
        if TIME.it==0 &&strcmpi(TRANSP.boundary_condition_start,'Constant')
            TRANSP.QS_start=TRANSP.QS(1);
        end
        if strcmpi(TRANSP.boundary_condition_start,'Closed') && ~COAST.cyclic 
            % TRANSP.QS(1)=0; % Condition handled in coastline_change
        elseif strcmpi(TRANSP.boundary_condition_start,'Neumann') && ~COAST.cyclic
            TRANSP.QS_start=TRANSP.QS(1);
        elseif strcmpi(TRANSP.boundary_condition_start,'Func') && ~COAST.cyclic
            BCraw=load(TRANSP.BCfile);warning off;
            BC=struct;
            BC.timenum=datenum([num2str(BCraw(:,1))],'yyyymmdd');
            BC.QS_strt = interpNANs(BCraw(:,2));
            BC.QSSi=interp1(BC.timenum,BC.QS_strt,TIME.tnow);
            TRANSP.QS_start=interpNANs(BC.QSSi);
        elseif strcmpi(TRANSP.boundary_condition_start,'CTAN') && ~COAST.cyclic
            TRANSP.QS(1)=TRANSP.QS_start;
        elseif strcmpi(TRANSP.boundary_condition_start,'Periodic') && ~COAST.cyclic
            TRANSP.QS_start=TRANSP.QS(end-1);
        end
    end
    %% end point (right)
    if ~isempty(GROYNE.x)
        [idgroyne]=find(GROYNE.idcoast(:,1)==COAST.i_mc);
    end

    if isempty(idgroyne)  % Do not apply any of these transport boundary conditions 
                        % in case of a groyne on this side
        if TIME.it==0 &&strcmpi(TRANSP.boundary_condition_end,'Constant')
            TRANSP.QS_end=TRANSP.QS(end);
        end
        if strcmpi(TRANSP.boundary_condition_end,'Closed') && ~COAST.cyclic 
            TRANSP.QS_end=0; % Condition handled in coastline_change
        elseif strcmpi(TRANSP.boundary_condition_end,'Neumann') && ~COAST.cyclic 
            TRANSP.QS_end=TRANSP.QS(end-1);
        elseif strcmpi(TRANSP.boundary_condition_end,'Function') &&  strcmpi(TRANSP.boundary_condition_start,'FUNC') && ~COAST.cyclic
            BC.QS_end=interpNANs(BCraw(:,3));
            BC.QSEi=interp1(BC.timenum,BC.QS_end,TIME.tnow);
            TRANSP.QS_end=interpNANs(BC.QSEi);
        elseif strcmpi(TRANSP.boundary_condition_end,'Function') &&  ~strcmpi(TRANSP.boundary_condition_start,'FUNC') && ~COAST.cyclic
            BCraw=load(TRANSP.BCfile);warning off;
            BC=struct;
            BC.timenum=datenum([num2str(BCraw(:,1))],'yyyymmdd');
            BC.QS_end=interpNANs(BCraw(:,3));
            BC.QSEi=interp1(BC.timenum,BC.QS_end,TIME.tnow);
            TRANSP.QS_end=interpNANs(BC.QSEi);
        elseif strcmpi(TRANSP.boundary_condition_end,'Constant') && ~COAST.cyclic 
            %TRANSP.QS_end=TRANSP.QS_end;
        elseif strcmpi(TRANSP.boundary_condition_end,'Periodic') && ~COAST.cyclic 
            TRANSP.QS_end=TRANSP.QS(2);
        end
    end

    % set nans to 0
    if find(isnan(TRANSP.QS))>0
        TRANSP.QS(find(isnan(TRANSP.QS)))=0;
    end

end

