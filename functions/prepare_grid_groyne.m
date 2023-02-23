function [COAST,STRUC,GROYNE,TRANSP]=prepare_grid_groyne(COAST,STRUC,TRANSP,yesplot)
% function [COAST,STRUC,GROYNE,TRANSP]=prepare_grid_groyne(COAST,STRUC,TRANSP,yesplot)
% 
% Initialize coastline grid taking into account groynes that intersect
% the coastline.
%
% Input: 
%   COAST
%        .x_mc      : coastline x-coordinates
%        .y_mc      : coastline y-coordinates
%        .ds0       : grid resolution (m)
%   STRUC
%        .x_hard    : hard structure elements. Groynes must be represented
%        .y_hard    : by polylines intersecting coastline twice; clockwise
%   yesplot         : switch for plotting
%
% Output:
%   COAST
%        .x_mc      : coastline x-coordinates
%        .y_mc      : coastline y-coordinates
%   STRUC
%        .s_hard    : cumulative distance along each groyne
%   GROYNE
%        .s         : along-groyne coordinates of both intersections for each groyne
%        .x         : x-coordinates for each groyne and intersection
%        .y         : y-coordinates for each groyne and intersection
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
    
    %fprintf('  Initialize groynes (grid) \n');   
    %% insert  points at the intersection between structure and coastline if exists  ( Added by M.Ghonim )
    eps=1e-1;
    GROYNE.s=[];
    GROYNE.x=[];
    GROYNE.y=[];
    GROYNE.n=0;
    GROYNE.QS=[];
    GROYNE.strucnum=[];
    GROYNE.idcoast=[];
                
    %% Find intersections with groynes, split coastline and store groyne intersections
    if length(STRUC.x_hard)>1
        n_st=length(find(isnan(STRUC.x_hard)))+1;                                    % Number of groynes
        ii=0;
        for ist=1:n_st
            [ xs,ys,~,~,~ ] = get_one_polygon( STRUC.x_hard,STRUC.y_hard,ist);
            ss=zeros(size(xs));
            for i=2:length(xs)
                ss(i)=ss(i-1)+hypot(xs(i)-xs(i-1),ys(i)-ys(i-1));
            end
            if ist==1
                STRUC.s_hard=[ss];
            else
                STRUC.s_hard=[STRUC.s_hard,nan,ss];
            end
 
            % find intersections of current structure with full coastline
            [str_int_x,str_int_y,indc,inds]=intersectfn(COAST.x_mc,COAST.y_mc,xs,ys);
            
            % correct the crossing points of the groyne with the coastline if it detects some twice (then use first and last one)
            [~,ids]=sort(indc);
            [~,idu]=unique(str_int_y.^2+str_int_x.^2);idu=idu(:)';
            if length(idu)==1 && length(indc)>=2
                indc=[indc(1),indc(2)+1];
                inds=[];
                str_int_x(2)=COAST.x_mc(indc(1));
                str_int_y(2)=COAST.y_mc(indc(2));
            elseif length(idu)>2
                ids=ids(idu);
            end
            if length(str_int_x)>2
                %disp('More than 2 crossing points of the coast and groyne are found! Only the first and last point are used.')
                ids=[ids(1),ids(end)];
            end
            indc=indc(ids)';
            inds=inds(ids)';
            str_int_x=str_int_x(ids);
            str_int_y=str_int_y(ids);
            
            % determine the x and y location where the groyne attaches to the coast (in GROYNE.x and GROYNE.y)
            % and store the coastal segment numbers the groyne is connected to (in GROYNE.idcoast)
            % then cut the coast in two pieces at both sides of the groyne. 
            if length(str_int_x) == 2
                % valid groyne, with two intersections
                ii=ii+1;
                GROYNE.n=GROYNE.n+1;
                GROYNE.strucnum(ii)=ist;
                GROYNE.s(ii,:)=interp1([1:length(ss)],ss,[min(inds),max(inds)]);
                GROYNE.x(ii,:)=str_int_x+diff(str_int_x)*[0.00001,-0.00001];
                GROYNE.y(ii,:)=str_int_y+diff(str_int_y)*[0.00001,-0.00001];
                GROYNE.xs=xs;
                GROYNE.ys=ys;
                GROYNE.QS=zeros(GROYNE.n,2);
                GROYNE.tipindx=zeros(GROYNE.n,2);
                GROYNE.phicoast(ii,1)=mod(360-atan2d(diff(str_int_y),diff(str_int_x)),360);
                dxgroyne=mean(xs(2:end))-(xs(1)+xs(end))/2;
                dygroyne=mean(ys(2:end))-(ys(1)+ys(end))/2;
                GROYNE.phigroyne(ii,1)=mod(90-atan2d(dygroyne,dxgroyne),360);
                
                % find segments index (i_mc) at both sides of the groyne split
                idnotnans=find(~isnan(COAST.x_mc(1:floor(indc(1)))));  % find the number of segments before         
                nrsegmentsbefore=length([idnotnans(diff(idnotnans)~=1)])+1; % count the number of segments 
                
                % add index of segments at both sides of the groyne
                GROYNE.idcoast(ii,:)=[nan,nan];
                GROYNE.idcoast(ii,:)=nrsegmentsbefore+[0,1];
                if GROYNE.idcoast
                if isempty(idnotnans)
                    GROYNE.idcoast(ii,1)=nan;
                end
                % check if there are segments after
                if sum(~isnan(COAST.x_mc(ceil(indc(2)):end)))==0 
                    GROYNE.idcoast(ii,2)=nan;
                end
                
                % make sure the segment-ids of earlier defined groyne-segment splits are adjusted (for segments after the current split)
                idelementsaftersplit=setdiff(find(GROYNE.idcoast(:,1)>=nrsegmentsbefore),ii);
                GROYNE.idcoast(idelementsaftersplit,:)=GROYNE.idcoast(idelementsaftersplit,:)+1;  
                
                % add nans to coastline
                COAST.x_mc = [COAST.x_mc(1:floor(indc(1))), GROYNE.x(ii,1),NaN,GROYNE.x(ii,2),COAST.x_mc(ceil(indc(2)):end)];
                COAST.y_mc = [COAST.y_mc(1:floor(indc(1))), GROYNE.y(ii,1),NaN,GROYNE.y(ii,2),COAST.y_mc(ceil(indc(2)):end)];
                               
                % add groyne number to array to make sure that the splitting nans of groynes can be traced back
                COAST.i_groyne(floor(indc(1))+1)=ist;  
                
            elseif length(str_int_x)<2
                disp('Less than 2 crossing points of the coast and groyne are found!');
                disp('Structure is assumed to be burried below the sand or acts as an offshore breakwater.');
            else
                disp(num2str(str_int_x))
            end  
        end
    end
    
    % determine number of coastline segments
    nans=find(isnan(COAST.x_mc));
    COAST.n_mc=length(nans)+1;
    
    if yesplot
       figure(11);
    end
end
