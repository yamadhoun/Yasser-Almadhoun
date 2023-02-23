function [x,y,xmax,ymax,xmin,ymin,nmax] = plot_insert_land_fill(x,y,ld,it,xmax,ymax,xmin,ymin,i_mc,n_mc,nmax)
% function [x,y,xmax,ymax,xmin,ymin,nmax] = plot_insert_land_fill(x,y,ld,it,xmax,ymax,xmin,ymin,i_mc,n_mc,nmax)
% insert land area behind shoreline , the area width ld (m)
% surrounding with black line
% work with fill_sections function
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
  
    if it==0  % this value should store
        xmin(i_mc)=min(x);
        ymin(i_mc)=min(y);
        xmax(i_mc)=max(x);
        ymax(i_mc)=max(y);
        nmax=n_mc;
    elseif i_mc>length(xmin)
        xmin(i_mc)=min(x);
        ymin(i_mc)=min(y);
        xmax(i_mc)=max(x);
        ymax(i_mc)=max(y);
        nmax=n_mc;
    end
    %  if n_mc > nmax
    %      nads=n_mc-nmax
    %      if i_mc > nmax
    %          for ina=1:nads
    %              xmin(i_mc)=min(x);
    %              ymin(i_mc)=min(y);
    %              xmax(i_mc)=max(x);
    %              ymax(i_mc)=max(y);
    %          end
    %      end 
    %  end
           
      
    o=atan2d((y(end)-y(1)),(x(end)-x(1)));
    dist=((y(end)-y(1)).^2+(x(end)-x(1)).^2).^0.5;
    distmax=max( max(((y-y(1)).^2+(x-x(1)).^2).^0.5), max(((y-y(end)).^2+(x-x(end)).^2).^0.5));
    ds0=hypot(x(2)-x(1),y(2)-y(1));
    
    % close the polygon
    if dist>0.2*distmax % may be check that dist>2*ds0 ?                %&& n_mc==1
        % close the polygon when it is quite open (using additional points of extents of the model)
        % otherwise close the polygon directly with begin/end points (when relatively close i.e. within 20% of max distance between points)
        if o >= -45 &&  o<=45   %% Shoreline horizontal & land on the right side (down)
             x=[xmin(i_mc),x,xmax(i_mc)];
             y=[min(ymin)-ld,y,min(ymin)-ld];
        elseif abs (o) >= 135  %% Shoreline horizontal & land on the left side (up)
             x=[xmax(i_mc),x,xmin(i_mc)];
             y=[max(ymax)+ld,y,max(ymax)+ld];
        elseif o > 45 &&  o<135    %% vertical Shoreline  & land on the right side 
            x=[max(xmax)+ld,x,max(xmax)+ld];
             y=[ymin(i_mc),y,ymax(i_mc)];
        elseif o < -45 &&  o>-135   %% vertical Shoreline  & land on the left side 
             x=[min(xmin)-ld,x,min(xmin)-ld];
             y=[ymax(i_mc),y,ymin(i_mc)];
        end
    end
    xb=[x(end-1),x(end),x(1),x(2)];
    yb=[y(end-1),y(end),y(1),y(2)];
       
    plot(xb,yb,'k','linewidth',3);
    hold on 
    
end
