function [DUNE]=prepare_dunes(S,COAST)
% function [DUNE]=prepare_dunes(S,COAST)
% 
% INPUT
%    S
%
% OUTPUT
%    DUNE
%       .x_dune
%       .y_dune
%       .x_Df
%       .y_Df
%       .Dfelev0
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

global intersectfn

    fprintf('  Prepare dunes \n');
    DUNE=struct;
    DUNE.used=S.dune;
    DUNE.x_dune=[];
    DUNE.y_dune=[];
    DUNE.x_Df=[];
    DUNE.y_Df=[];
    DUNE.Dfelev0=[]; 
    DUNE.Dfelevation=[];
    DUNE.x_dune0=[];
    DUNE.y_dune0=[];  
    DUNE.Dfelev=[];
    
    %% PREPARE DUNE PROPERTIES
    if DUNE.used
        if ~isempty(S.LDBdune)&&~isfield(S,'DUNE.x_dune')
                xDUNE.y_dune=load(S.LDBdune);
            if isempty(S.XYoffset)
                S.XYoffset = [floor(min(xDUNE.y_dune(:,1))/1000)*1000 , floor(min(xDUNE.y_dune(:,2))/1000)*1000];
            end
            DUNE.x_dune=xDUNE.y_dune(:,1)' - S.XYoffset(1);   
            DUNE.y_dune=xDUNE.y_dune(:,2)' - S.XYoffset(2);   

        elseif ~isfield(S,'DUNE.x_dune')
            figure(11);
            axis equal;
            xl=xlim;yl=ylim;
            htxt=text(xl(1)+0.02*diff(xl),yl(2)-0.01*diff(yl),'Add Dune foot (LMB); Next segment (RMB); Exit (q)');set(htxt,'HorizontalAlignment','Left','VerticalAlignment','Top','FontWeight','Bold','FontAngle','Italic','Color',[0.6 0.1 0.1]);
            [DUNE.x_dune,DUNE.y_dune]=select_multi_polygon('g');
            set(htxt,'Visible','off');

        else
            DUNE.x_dune=S.DUNE.x_dune;
            DUNE.y_dune=S.DUNE.y_dune;

        end
        if~isempty(S.Dfelevation)
            xDUNE.y_Df=load(S.Dfelevation);
            DUNE.x_Df=xDUNE.y_Df(:,1)';
            DUNE.y_Df=xDUNE.y_Df(:,2)';
            DUNE.Dfelev0=xDUNE.y_Df(:,3)';
        else
            DUNE.x_Df=[];
            DUNE.y_Df=[];
            DUNE.Dfelev0=S.Dfelev;
        end
        DUNE.Dfelevation=S.Dfelevation;

    %% PREPARE DUNE GRID
        phid=zeros(size(COAST.x_mc));
        len=hypot(max(COAST.x_mc)-min(COAST.x_mc),max(COAST.y_mc)-min(COAST.y_mc));
        nansd=find(isnan(DUNE.x_dune));
        n_dune=length(nansd)+1;
        for i_dune=1:n_dune
            [ xd0,yd0,~,id1,id2 ] = get_one_polygon( DUNE.x_dune,DUNE.y_dune,i_dune );
            for i=1:length(COAST.x_mc)-1
                phic(i)=360-atan2d(COAST.y_mc(i+1)-COAST.y_mc(i),COAST.x_mc(i+1)-COAST.x_mc(i));
            end
            for j=2:length(COAST.x_mc)-1
                phid(j)=0.5*(phic(j)+phic(j-1));
            end
            phid(1)=phic(1);
            phid(end)=phic(end);
            for ii=1:length(COAST.x_mc)
                Xs=[COAST.x_mc(ii)-len*sind(phid(ii)),COAST.x_mc(ii)+len*sind(phid(ii))]; 
                Ys=[COAST.y_mc(ii)-len*cosd(phid(ii)),COAST.y_mc(ii)+len*cosd(phid(ii))];
                [xx1,yy1]=intersectfn(Xs,Ys,xd0,yd0);
                if ~isempty(xx1) && length(xx1)==2
                    if hypot(COAST.x_mc(ii)-xx1(1),COAST.y_mc(ii)-yy1(1))<hypot(COAST.x_mc(ii)-xx1(2),COAST.COAST.y_mc(ii)-yy1(2))
                        xd(ii)=xx1(1);
                        yd(ii)=yy1(1);
                    else 
                        xd(ii)=xx1(2);
                        yd(ii)=yy1(2);
                    end
                elseif ~isempty(P1)
                    xd(ii)=xx1(1);
                    yd(ii)=yy1(1);
                else
                     xd(ii)=-1;
                     yd(ii)=-1;
                end
            end
            ndx=xd==-1;
            ndy=yd==-1;
            xd(ndx)=[];
            yd(ndy)=[];     
            %% insert xd,yd back into DUNE.x_dune,DUNE.y_dune
            DUNE.x_dune=[DUNE.x_dune(1:id1-1),xd,DUNE.x_dune(id2+1:end)];
            DUNE.y_dune=[DUNE.y_dune(1:id1-1),yd,DUNE.y_dune(id2+1:end)];
            xy_dune(1,:)=DUNE.x_dune;
            xy_dune(2,:)=DUNE.y_dune;
            for aa=2:length(xy_dune)
                if all(xy_dune(:,aa)==xy_dune(:,(aa-1)))            % remove any duplicated full columns if exist
                    xy_dune(:,aa)=[];
                end
            end
            DUNE.x_dune=xy_dune(1,:);
            DUNE.y_dune=xy_dune(2,:);
            DUNE.x_dune0=DUNE.x_dune;
            DUNE.y_dune0=DUNE.y_dune;
            xy_dune=[];
            % interpolate for dune foot height at DUNE.x_dune & DUNE.y_dune locations
            if~isempty(DUNE.Dfelevation)
                s0=zeros(size(x_Df));
                s=zeros(size(DUNE.x_dune));
                DUNE.Dfelev=Dfelev0;
                for ii=2:length(x_Df)
                    s0(ii)=s0(ii-1)+hypot(x_Df(ii)-x_Df(ii-1),y_Df(ii)-y_Df(ii-1));
                end
                for i=2:length(DUNE.x_dune)
                    s(i)=s(i-1)+hypot(DUNE.x_dune(i)-DUNE.x_dune(i-1),DUNE.y_dune(i)-DUNE.y_dune(i-1));
                end
                DUNE.Dfelev=interp1(s0,DUNE.Dfelevation,s);
            else
                DUNE.Dfelev=Dfelev0.*ones(1,length(DUNE.x_dune));
            end
        end
    end
end