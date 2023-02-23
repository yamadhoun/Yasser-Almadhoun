function [xl,yl,phirf]= dune_refline(S,x_dune,y_dune) 
% function [xl,yl,phirf]= dune_refline(S,x_dune,y_dune) 
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
    xl=zeros(size(x_dune));
    yl=zeros(size(x_dune));
    nansd=find(isnan(x_dune));
    n_dune=length(nansd)+1;
    for i_dune=1:n_dune
        [ xd,yd,~,id1,id2 ] = get_one_polygon( x_dune,y_dune,i_dune );
    %     cyclicd = hypot(x_dune(end)-x_dune(1),y_dune(end)-y_dune(1))<S.ds0;
    %     if cyclicd
    %         nend=nd;
    %     else
        nd=length(xd)-1;
        nend=nd+1;
        phirf0=zeros(size(xd));
    %     end
        for i=1:nend
    %         if cyclicd
    %             im1=mod2(i-1,nd);
    %             ip1=mod2(i+1,nd);
    %         else
                im1=max(i-1,1);
                ip1=min(i+1,nd+1);
                % computed angles range from 180 to 540°
                phirf0(i)=360-atan2d(yd(ip1)-yd(im1),xd(ip1)-xd(im1));
    %         end
        end
        phirf(id1:id2)=phirf0;
    end
    for ii=1:length(x_dune)
        Xl=[x_dune(ii),x_dune(ii)-500*sind(phirf(ii))]; 
        Yl=[y_dune(ii),y_dune(ii)-500*cosd(phirf(ii))];
        xl(ii)=Xl(end);
        yl(ii)=Yl(end);
    end
else
    xl=[];
    yl=[];
    phirf=[];
end
end