function [COAST]=prepare_coastline(S)
% function [x_mc,y_mc,x_mc0,y_mc0,S]=prepare_coastline(S)
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

    fprintf('  Prepare coastline \n');
    COAST=struct;
    COAST.LDBcoastline=S.LDBcoastline;
    COAST.h0=S.d;
    COAST.ds0=S.ds0;
    COAST.xlimits=S.xlimits;
    COAST.ylimits=S.ylimits;
    COAST.XYoffset=S.XYoffset;
    COAST.twopoints=S.twopoints;
    COAST.smoothfac=S.smoothfac;
    COAST.PHIf0=S.phif;
    COAST.tanbeta=S.tanbeta;
    
    if ~isempty(S.LDBcoastline) && ~isfield(S,'x_mc')
        if ischar(S.LDBcoastline)
            xy_mc=load(S.LDBcoastline);
        else
            xy_mc=S.LDBcoastline;
        end
        if isempty(S.XYoffset)
            COAST.XYoffset = [floor(min(xy_mc(:,1))/1000)*1000 , floor(min(xy_mc(:,2))/1000)*1000];
        end
        COAST.x_mc=xy_mc(:,1)' - S.XYoffset(1);   %COAST.x_mc=xy_mc(end:-1:1,1)';   % SHIFT COASLTINE
        COAST.y_mc=xy_mc(:,2)' - S.XYoffset(2);   %COAST.y_mc=xy_mc(end:-1:1,2)';   % SHIFT COASLTINE
        COAST.x_mc0=COAST.x_mc;
        COAST.y_mc0=COAST.y_mc;
        
    elseif ~isfield(S,'x_mc')
        figure(11);clf;
        %plot_figureproperties(gcf,800,950,32);
        %plot_figureproperties(gcf,920,1120,32,80,0);
        %plot_figureproperties(gcf,1855,1120,32,65,0);
        xl=S.xlimits;yl=S.ylimits;
        plot([S.xlimits(1) S.xlimits(2) S.xlimits(2) S.xlimits(1) S.xlimits(1)], ...
            [S.ylimits(1) S.ylimits(1) S.ylimits(2) S.ylimits(2) S.ylimits(1)],'k:');
        axis equal;
        %xl=xlim;yl=ylim;
        xlabel('Easting [m]');
        ylabel('Northing [m]');
        htxt=text(xl(1)+0.02*diff(xl),yl(2)-0.01*diff(yl),'Add coastline (LMB); Next segment (RMB); Exit (q)');set(htxt,'HorizontalAlignment','Left','VerticalAlignment','Top','FontWeight','Bold','FontAngle','Italic','Color',[0.1 0.1 0.6]);
        [COAST.x_mc,COAST.y_mc]=select_multi_polygon('r');
        set(htxt,'Visible','off');
        COAST.x_mc0=COAST.x_mc;
        COAST.y_mc0=COAST.y_mc;

    else
        COAST.x_mc=S.x_mc;
        COAST.y_mc=S.y_mc;
        try
            COAST.x_mc0=S.x_mc0;
            COAST.y_mc0=S.y_mc0;
        catch 
            COAST.x_mc0=S.x_mc;
            COAST.y_mc0=S.y_mc;
        end
        figure(11);clf;
    end
    
    % tidy up the coastline without nans at the start and end, and remove double nan's
    idnotnan=find(~isnan(COAST.x_mc));
    COAST.x_mc=COAST.x_mc(idnotnan(1):idnotnan(end));
    COAST.y_mc=COAST.y_mc(idnotnan(1):idnotnan(end));
    idnan=find(isnan(COAST.x_mc));
    iduse=setdiff([1:length(COAST.x_mc)],idnan(diff(idnan)==1));
    COAST.x_mc=COAST.x_mc(iduse);
    COAST.y_mc=COAST.y_mc(iduse);
    % determine number of coastline segments
    nans=find(isnan(COAST.x_mc));
    COAST.n_mc=length(nans)+1;
    
end 
