function [ shadowS ] = find_shadows_mc(xq,yq,x_mc,y_mc,PHI,PHItdp,PHIbr)
%
% Computes the indices of locations in the shadow zone 
% based on coastline shape, location of structures and wave incidence angle.
%
% INPUT:
%         x         : x-coordinate of coastline (only current section)
%         y         : y-coordinate of coastline (only current section)
%         x_mc      : x-coordinate of coastline (all sections)
%         y_mc      : y-coordinate of coastline (all sections)
%         PHI       : offshore wave incidence angle ([1] or [1xN] in degrees North)
%         PHItdp    : wave incidence angle at depth of closure ([1] or [1xN] in degrees North, using offshore wave if CERC is used)
%         PHIbr     : wave incidence angle at breaking point ([1] or [1xN] in degrees North, using offshore wave if CERC is used)
%         hard      : switch for hard structures
%
% OUTPUT:
%         xS        : x-coordinate of QS-points
%         yS        : y-coordinate of QS-points
%         shadowS   : Index of cells which are in the shadow zone
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
    
    nq=length(PHItdp);   %length(x)-1;
    if length(PHI)==1
        PHI=repmat(PHI,[1,nq]);
    end
    lenSURF=300;
    lenBREAKER=100;
    
    if nq==0
        shadowS=[];
    else
        len=5*hypot(max(xq)-min(xq),max(yq)-min(yq));
        shadowS=false(1,nq);
        for i=1:nq
            % using a line width a small bend at the point of breaking (using PHIbr) and nearshore depth of closure (using PHItdp)
            xw=[xq(i)+1*sind(PHIbr(i)),...
                xq(i)+lenBREAKER*sind(PHIbr(i)),...
                xq(i)+lenBREAKER*sind(PHIbr(i))+lenSURF*sind(PHItdp(i)),...
                xq(i)+lenBREAKER*sind(PHIbr(i))+lenSURF*sind(PHItdp(i))+len*sind(PHI(i))];
            yw=[yq(i)+1*cosd(PHIbr(i)),...
                yq(i)+lenBREAKER*cosd(PHIbr(i)),...
                yq(i)+lenBREAKER*cosd(PHIbr(i))+lenSURF*cosd(PHItdp(i)),...
                yq(i)+lenBREAKER*cosd(PHIbr(i))+lenSURF*cosd(PHItdp(i))+len*cosd(PHI(i))];
            
            [xx1,~,~,~,~,~,~,~]=get_intersections(x_mc,y_mc,xw,yw);
            
            shadowS(i)=~isempty(xx1);
            if 0
                figure(10);clf
                % plot(xq,yq,x_mc,y_mc,xw,yw,'.-',xx1,yy1,'ro','linewidth',2);
                plot(x,y,'b+');hold on;
                plot(x_mc,y_mc,'r.-');
                plot(xw,yw,'g.-');
                plot(P1(1,:),P1(2,:),'ro','linewidth',1);
                xlim([min(x),max(x)]);
                ylim([min(y),max(y)]);
                axis equal
                drawnow
                %%pause
            end
        end
    end
end
