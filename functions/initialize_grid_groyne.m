function [COAST,STRUC,GROYNE,TRANSP]=initialize_grid_groyne(COAST,STRUC,TRANSP,yesplot)
% function [COAST,STRUC,GROYNE]=initialize_grid(COAST,STRUC,yesplot)
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
    eps=1e-6;
    GROYNE.s=[];
    GROYNE.x=[];
    GROYNE.y=[];
    GROYNE.n=0;
    GROYNE.QS=[];
    
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
 
            % Find intersections of current structure with full coastline
            [str_int_x,str_int_y, indc, inds] = intersectfn(COAST.x_mc,COAST.y_mc,xs,ys);
            if length(str_int_x)>2
                IDunique = [true,hypot(diff(str_int_x),diff(str_int_y))>eps];
                str_int_x = str_int_x(IDunique);
                str_int_y = str_int_y(IDunique);
                indc = indc(IDunique);
                inds = inds(IDunique);
            end
            if length(str_int_x) == 2
                % valid groyne, with two intersections
                ii=ii+1;
                GROYNE.strucnum(ii)=ist;
                GROYNE.s(ii,:)=interp1([1:length(ss)],ss,inds);
                GROYNE.x(ii,:)=str_int_x;
                GROYNE.y(ii,:)=str_int_y;
                COAST.x_mc = [COAST.x_mc(1:floor(indc(1))), GROYNE.x(ii,1),NaN,GROYNE.x(ii,2),COAST.x_mc(ceil(indc(2)):end)];
                COAST.y_mc = [COAST.y_mc(1:floor(indc(1))), GROYNE.y(ii,1),NaN,GROYNE.y(ii,2),COAST.y_mc(ceil(indc(2)):end)];
            end
                       
        end
    end
    % if (find(isnan(COAST.x_mc))~= COAST.n_mc)
    %     error('Initialize coastline:: Number of coastline sections is not equal to Num of structures +1 .');
    % end
    nans=find(isnan(COAST.x_mc));
    COAST.n_mc=length(nans)+1;

    if yesplot
       figure(11);
    end

    if 1
    %% Generate grid with uniform distances at approx. grid size S.COAST.ds0
    for i_mc=1:COAST.n_mc
        [ x,y,COAST.n_mc,i1,i2 ] = get_one_polygon( COAST.x_mc,COAST.y_mc,i_mc );   
        if yesplot
           plot(x,y); hold on
        end
        % remove double points
        IDunique = [true,hypot(diff(x),diff(y))>eps];
        x0=x(IDunique);
        y0=y(IDunique);
        % compute distance along line
        s0=zeros(size(x0));
        for i=2:length(x0)
            s0(i)=s0(i-1)+hypot(x0(i)-x0(i-1),y0(i)-y0(i-1));
        end
        ns=ceil(s0(end)/COAST.ds0);
        ds=s0(end)/ns;
        s=[0:ds:s0(end)];
        % interpolate x-values
        if length(s)>1
            try
            x=interp1(s0,x0,s);
            y=interp1(s0,y0,s);
            catch
                x
            end
            %% insert x,y back into COAST.x_mc,COAST.y_mc
            COAST.x_mc=[COAST.x_mc(1:i1-1),x,COAST.x_mc(i2+1:end)];
            COAST.y_mc=[COAST.y_mc(1:i1-1),y,COAST.y_mc(i2+1:end)];
        else
            COAST.n_mc=1;
        end
    end
    COAST.n_mc=length(find(isnan(COAST.x_mc)))+1;
    if yesplot
       plot(COAST.x_mc,COAST.y_mc,COAST.x_mc,COAST.y_mc,'.',STRUC.x_hard,STRUC.y_hard,'k')
       axis equal
    end
    end
    
    GROYNE.n=size(GROYNE.x,1);
    GROYNE.QS=zeros(GROYNE.n,2);
    
end
