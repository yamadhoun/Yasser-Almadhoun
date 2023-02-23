function [FORMAT] = initialize_plot(S,COAST)
% function [FORMAT.plot_time,FORMAT.tsl,FORMAT.tplot,FORMAT.xp_mc,FORMAT.yp_mc] = initialize_plot_variables(S)
% function [xl,yl,phirf,vi,vii,FORMAT.xmax,FORMAT.ymax,FORMAT.xmin,FORMAT.ymin,FORMAT.nmax,iwtw,FORMAT.plot_time,FORMAT.tsl,FORMAT.tplot,FORMAT.xp_mc,FORMAT.yp_mc,CLplot,CLplot2,iint,BWplot,BWplot2,innt,ds_cl,qwave,qwind,time,step,int,bermW,x_trans,y_trans,n_trans] = initialize_plot_variables(S,x_mc,y_mc,x_hard,y_hard,x_dune,y_dune,timenum0)
%
% UNTITLED5 Summary of this function goes here
% Detailed explanation goes here
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

    fprintf('  Initialize plot \n');
    
    % Format figures
    if S.plotvisible==0
        S.plotvisible='off';
    else
        S.plotvisible='on';
    end
    FORMAT.mainfighandle = figure(11);clf;
    set(FORMAT.mainfighandle,'Visible',S.plotvisible);
    plot_figureproperties(FORMAT.mainfighandle,800,600,32);
    FORMAT.xmax=0;FORMAT.ymax=0;FORMAT.xmin=0;FORMAT.ymin=0;FORMAT.nmax=0;
    FORMAT.plotvisible=S.plotvisible;
    
    % Set x and y limits
    if isempty(S.xlimits)
        S.xlimits = [min(COAST.x_mc),max(COAST.x_mc)];
    else
        S.xlimits = S.xlimits - S.XYoffset(1);
    end
    if isempty(S.ylimits)
        S.ylimits = [min(COAST.y_mc),max(COAST.y_mc)];
    else
        S.ylimits = S.ylimits - S.XYoffset(2);
    end
    S.xlimits=sort(S.xlimits);
    S.ylimits=sort(S.ylimits);
    xlim(S.xlimits);ylim(S.ylimits);
    FORMAT.xlimits=S.xlimits;
    FORMAT.ylimits=S.ylimits;
    FORMAT.XYoffset=S.XYoffset;

    if S.yesplot
        figure(11);clf;
        plot([S.xlimits(1) S.xlimits(2) S.xlimits(2) S.xlimits(1) S.xlimits(1)], ...
            [S.ylimits(1) S.ylimits(1) S.ylimits(2) S.ylimits(2) S.ylimits(1)],'k:');
        axis equal;
        hold on;
        xlabel('Easting [m]');
        ylabel('Northing [m]');
        plot(COAST.x_mc0,COAST.y_mc0,'k','linewidth',2);
        axis equal
        if ~isempty(S.xlimits)
            xlim(S.xlimits);
        end
        if ~isempty(S.ylimits)
            ylim(S.ylimits);
        end
    end

    FORMAT.xp_mc={};
    FORMAT.yp_mc={};

    % Other plot formatting
    if ~isempty(S.SLplot)  % For extracting specific shorelines
        FORMAT.plot_time(:)=datenum(S.SLplot(:,1),'yyyy-mm-dd');
        FORMAT.tsl=1;
        FORMAT.tplot=FORMAT.plot_time(FORMAT.tsl);
        FORMAT.plot_time(end+1)=0;
    else
        FORMAT.tplot=[];
        FORMAT.tsl=[];
        FORMAT.plot_time=[];
    end
    
    FORMAT.plotinterval=S.plotinterval;
    FORMAT.usefill=S.usefill;
    FORMAT.LDBplot=S.LDBplot;
    FORMAT.llocation=S.llocation;
    FORMAT.SLplot=S.SLplot;
    FORMAT.video=S.video;
    FORMAT.outputdir=S.outputdir;
    FORMAT.fignryear=S.fignryear;
    FORMAT.ld=S.ld;
    FORMAT.XYwave=S.XYwave;
    FORMAT.plotQS=S.plotQS;
    FORMAT.print_fig=S.print_fig;
    FORMAT.fastplot=S.fastplot;
end
 
% %  Create reference line for cross-shore plotting
% [xl,yl,phirf]= dune_refline(S,x_dune,y_dune);
% 
% % for making video
% vi = struct('cdata', cell(1,1), 'colormap', cell(1,1));
% vii=0;
% 
% %for plotting fill
% FORMAT.xmax=0;FORMAT.ymax=0;FORMAT.xmin=0;FORMAT.ymin=0;FORMAT.nmax=0;
% iwtw=0;  % for using interpolated wave table (warning)
% 
% % For extracting specific shorelines
% if ~isempty(S.SLplot)
%     FORMAT.plot_time(:)=datenum(S.SLplot(:,1),'yyyy-mm-dd'); % HH:MM:SS
%     FORMAT.tsl=1;
%     FORMAT.tplot=FORMAT.plot_time(FORMAT.tsl);
%     FORMAT.plot_time(end+1)=0;
% else
%     FORMAT.tplot=[];
%     FORMAT.tsl=[];
%     FORMAT.plot_time=[];
% end
% FORMAT.xp_mc{:,:}={};
% FORMAT.FORMAT.yp_mc{:,:}={};
% if S.dune && S.qplot                                                         % for extracting wave and wind transport rates at a certain transect(s)
%     innt=1;
% else
%     innt=[];
% end
% if S.dune && S.bermw_plot 
%     BWplot=timenum0;                                                       % for extracting beach berm width at a certain transect(s)
%     BWplot2=timenum0;
%     int=1;
% else
%     BWplot=[];
%     BWplot2=[];
%     int=[];
% end
% if S.CLplot
%     CLplot=timenum0;                                                       % for extracting coastline position relative to the initial coastline at a certain transect(s)
%     CLplot2=timenum0;
%     iint=1;
% else
%     CLplot=[];
%     CLplot2=[];    
%     iint=1;
% end
% ds_cl=[];
% qwave=[];
% qwind=[];
% time=[];
% step=[];
% bermW=[];
% 
% %% prepare transections from cross-shore plotting
% if S.CLplot ||S.bermw_plot ||S.qplot 
%     [x_trans,y_trans,n_trans]=prepare_transects(S,x_mc,y_mc,x_hard,y_hard,x_dune,y_dune,xl,yl);
% else
%     x_trans=[];
%     y_trans=[];
%     n_trans=[];
% end
% end
% 
