function [x_dune,y_dune,x_dune0,y_dune0,Dfelev]=initialize_grid_dune(S,x_mc,y_mc,x_dune,y_dune,x_Df,y_Df,Dfelev0)
% function [x_dune,y_dune,x_dune0,y_dune0,Dfelev]=initialize_grid_dune(S,x_mc,y_mc,x_dune,y_dune,x_Df,y_Df,Dfelev0)
% 
% intialize dune foots as intersect with coasline normal
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

    if S.dune
        phid=zeros(size(x_mc));
        len=hypot(max(x_mc)-min(x_mc),max(y_mc)-min(y_mc));
        nansd=find(isnan(x_dune));
        n_dune=length(nansd)+1;
        for i_dune=1:n_dune
            [ xd0,yd0,~,id1,id2 ] = get_one_polygon( x_dune,y_dune,i_dune );
            for i=1:length(x_mc)-1
                phic(i)=360-atan2d(y_mc(i+1)-y_mc(i),x_mc(i+1)-x_mc(i));
            end
            for j=2:length(x_mc)-1
                phid(j)=0.5*(phic(j)+phic(j-1));
            end
            phid(1)=phic(1);
            phid(end)=phic(end);
            for ii=1:length(x_mc)
                Xs=[x_mc(ii)-len*sind(phid(ii)),x_mc(ii)+len*sind(phid(ii))]; 
                Ys=[y_mc(ii)-len*cosd(phid(ii)),y_mc(ii)+len*cosd(phid(ii))];
                [xx1,yy1]=intersectfn(Xs,Ys,xd0,yd0);
                if ~isempty(xx1) && length(xx1)==2
                    if hypot(x_mc(ii)-xx1(1),y_mc(ii)-yy1(1))<hypot(x_mc(ii)-xx1(2),y_mc(ii)-yy1(2))
                        xd(ii)=xx1(1);
                        yd(ii)=yy1(1);
                    else 
                        xd(ii)=xx1(2);
                        yd(ii)=yy1(2);
                    end
                elseif ~isempty(xx1)
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
            %% insert xd,yd back into x_dune,y_dune
            x_dune=[x_dune(1:id1-1),xd,x_dune(id2+1:end)];
            y_dune=[y_dune(1:id1-1),yd,y_dune(id2+1:end)];
            xy_dune(1,:)=x_dune;
            xy_dune(2,:)=y_dune;
            for aa=2:length(xy_dune)
                if all(xy_dune(:,aa)==xy_dune(:,(aa-1)))            % remove any duplicated full columns if exist
                    xy_dune(:,aa)=[];
                    if 1
                        break
                    end
                end
            end
            x_dune=xy_dune(1,:);
            y_dune=xy_dune(2,:);
            x_dune0=x_dune;
            y_dune0=y_dune;
            xy_dune=[];
            % interpolate for dune foot height at x_dune & y_dune locations
            if~isempty(S.Dfelevation)
                s0=zeros(size(x_Df));
                s=zeros(size(x_dune));
                Dfelev=Dfelev0;
                for ii=2:length(x_Df)
                    s0(ii)=s0(ii-1)+hypot(x_Df(ii)-x_Df(ii-1),y_Df(ii)-y_Df(ii-1));
                end
                for i=2:length(x_dune)
                    s(i)=s(i-1)+hypot(x_dune(i)-x_dune(i-1),y_dune(i)-y_dune(i-1));
                end
                Dfelev=interp1(s0,Dfelev,s);
            else
                Dfelev= Dfelev0.* ones(1,length(x_dune));
            end
        end
    else
        x_dune=[];
        y_dune=[];
        x_dune0=[];
        y_dune0=[];  
        Dfelev=[];
    end
end
