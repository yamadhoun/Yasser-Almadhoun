function [xmax,ymax,xmin,ymin,nmax ] = plot_fill_sections(x_mc,y_mc,ld,it,xmax,ymax,xmin,ymin,nmax,ds0,usefill)
% function [xmax,ymax,xmin,ymin,nmax ] = plot_fill_sections(x_mc,y_mc,ld,it,xmax,ymax,xmin,ymin,nmax,ds0,usefill)
%
% UNTITLED Summary of this function goes here
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

    if nargin<3
        ld=[];
        it=[];
        nmax=[];
    end
    if nargin<6
        xmax=[];
        ymax=[];
        xmin=[];
        ymin=[];
    end
    if nargin<11
        usefill=1;
    end
    
    nans=find(isnan(x_mc));
    n_mc=length(nans)+1;
    
    for i_mc=1:n_mc
        nans=find(isnan(x_mc));
        if isempty(nans)
            i1=1; 
            i2=length(x_mc);
        else
            if i_mc==1
                i1=1;
                i2=nans(i_mc)-1;
            elseif i_mc==n_mc;
                i1=nans(i_mc-1)+1;
                i2=length(x_mc);
            else
                i1=nans(i_mc-1)+1;
                i2=nans(i_mc)-1;
            end
        end
        x=x_mc(i1:i2);
        y=y_mc(i1:i2);
        x0{i_mc}=x;
        y0{i_mc}=y;
    %     phi=[];
    %     for i=1:length(x)-1
    %         phi(i)=atan2d((y(i+1)-y(i)),(x(i+1)-x(i)));
    %     end
    % 
    %     % MAKE SURE TO FILL LAND AND LAKES DIFFERENTLY
    %     dphi=[];
    %     for i=1:length(x)-2
    %         dphi(i)=phi(i+1)-phi(i);
    %         if dphi(i)>180
    %             dphi(i)=dphi(i)-360;
    %         end
    %         if dphi(i)<-180
    %             dphi(i)=dphi(i)+360;
    %         end
    %     end
    %     %right=sum(dphi)<(1.8*180);                                              % theoretically a lake should have : sum(dphi)==360 (but in practice not exactly)
    %     right(i_mc)=sum(dphi)<0;                                              % theoretically a lake should have : sum(dphi)==360 (but in practice not exactly)
        right(i_mc)=1;    
        if ~isempty(x)
            if hypot(x(end)-x(1),y(end)-y(1))>ds0
               % for line element coast -> use land
               right(i_mc)=1;
            else
               % for closed polygons use clockpoly
               right(i_mc)=min(get_clockpoly(x,y)+1,1);
            end
        end
    end
    IDlake=find(right==0);
    
    for i_mc=[setdiff([1:n_mc],IDlake),IDlake]
        x=x0{i_mc};
        y=y0{i_mc};
        notclosed=1;    % <- notclosed option recently added for Ameland (check?)
        if isempty(x)
            notclosed=0;
        elseif (x(1)-x(end)).^2+(y(1)-y(end)).^2 == 0
            notclosed=0;
        end
        
        if 0
        % CONSTRUCT A LAND POLYGON WHICH DOES NOT CROSS THE COASTLINE
        if (right(i_mc) || i_mc==1) && notclosed  % <- notclosed option recently added for Ameland (check?)
            dx=diff(x);dx=[dx(1),(dx(1:end-1)+dx(2:end))/2,dx(end)];
            dy=diff(y);dy=[dy(1),(dy(1:end-1)+dy(2:end))/2,dy(end)];
            l2=(dx.^2+dy.^2).^0.5;
            xclose = [x(end)+dy(end)/l2(end),x(1)+dy(1)/l2(1)];
            yclose = [y(end)-dx(end)/l2(end),y(1)-dx(1)/l2(1)];
            % plot(x_mc,y_mc,dx,dy,xclose,yclose,'*')
    
            % look for points left and right of line
            RC = diff(yclose)/diff(xclose);
            if diff(xclose)==0
               RC = diff(yclose)/1e-5;
            end
            B = yclose(1)-RC*xclose(1);
            ID = find(y<RC*x+B);                                       
            %ID = find(y>RC*x+B);
            %[xcr,ycr,x1n,y1n,x2n,y2n,ID]=findCrossings2(x,y,xclose,yclose);       % <- more sophisticated approach
            l3=100;
            xclose = [xclose(1),fliplr(x(ID)+l3*dy(ID)./l2(ID)),xclose(2)];
            yclose = [yclose(1),fliplr(y(ID)-l3*dx(ID)./l2(ID)),yclose(2)];
            x = [x,xclose];
            y = [y,yclose];
        end
        end
    
        if usefill==1
            % For straight line add land fill behind shoreline
            if notclosed && ~isempty(ld)
                [x,y,xmax,ymax,xmin,ymin,nmax] = plot_insert_land_fill(x,y,ld,it,xmax,ymax,xmin,ymin,i_mc,n_mc,nmax);
                %[x,y,nmax] = plot_insert_land_fill(x,y,ld,it,i_mc,n_mc,nmax);
            end
            % fill points
            if right(i_mc)
                fill(x,y,[.9 .9 0.4]);
            else
                fill(x,y,[.9 .9 1]);
            end
        else
            plot(x,y,'r','linewidth',2)
        end
        hold on
        
    end
    hold off

end
