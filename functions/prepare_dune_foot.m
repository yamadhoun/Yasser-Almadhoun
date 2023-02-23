function [x_dune,y_dune,x_Df,y_Df,Dfelev0]=prepare_dune_foot(S)
% function [x_dune,y_dune,x_Df,y_Df,Dfelev0]=prepare_dune_foot(S)
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

    if S.dune
        if ~isempty(S.LDBdune)&&~isfield(S,'x_dune')
                xy_dune=load(S.LDBdune);
            if isempty(S.XYoffset)
                S.XYoffset = [floor(min(xy_dune(:,1))/1000)*1000 , floor(min(xy_dune(:,2))/1000)*1000];
            end
            x_dune=xy_dune(:,1)' - S.XYoffset(1);   
            y_dune=xy_dune(:,2)' - S.XYoffset(2);   
    
        elseif ~isfield(S,'x_dune')
            figure(11);
            axis equal;
            xl=xlim;yl=ylim;
            htxt=text(xl(1)+0.02*diff(xl),yl(2)-0.01*diff(yl),'Add Dune foot (LMB); Next segment (RMB); Exit (q)');set(htxt,'HorizontalAlignment','Left','VerticalAlignment','Top','FontWeight','Bold','FontAngle','Italic','Color',[0.6 0.1 0.1]);
            [x_dune,y_dune]=select_multi_polygon('g');
            set(htxt,'Visible','off');
    
        else
            x_dune=S.x_dune;
            y_dune=S.y_dune;
    
        end
        if~isempty(S.Dfelevation)
            xy_Df=load(S.Dfelevation);
            x_Df=xy_Df(:,1)';
            y_Df=xy_Df(:,2)';
            Dfelev0=xy_Df(:,3)';
        else
            x_Df=[];
            y_Df=[];
            Dfelev0=S.Dfelev;
        end
    else
        x_dune=[];
        y_dune=[];
        x_Df=[];
        y_Df=[];
        Dfelev0=[]; 
    end
end
