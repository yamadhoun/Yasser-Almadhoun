function [x_dune,y_dune,S]=dune_foot_change(S,x_dune,y_dune,qss,qww,str_indx,thetaS)
% function [x_dune,y_dune,S]=dune_foot_change(S,x_dune,y_dune,qss,qww,str_indx,thetaS)
% 
% Estimate dune foot change as a result of wave action and aelian transport
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if S.dune 
    nansd=find(isnan(x_dune));
    n_dune=length(nansd)+1;
    for i_dune=1:n_dune
        [ xd,yd,~,id1,id2 ] = get_one_polygon( x_dune,y_dune,i_dune );
        qww_new=qww(id1:id2);
        qss_new=qss(id1:id2);
        dSdn=zeros(size(xd));
        dndt=zeros(size(xd));
        dxd=zeros(size(xd));
        dyd=zeros(size(xd));  
        nd=length(xd);
        cyclicd = hypot(xd(end)-xd(1),yd(end)-yd(1))<S.ds0;
   for i=1:nd
        dSdn(i)=(qww_new(i)-qss_new(i));
        dndt(i)=dSdn(i)/S.dn;
   end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dune foot change
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if cyclicd
        nend=nd-1;
    else
        nend=nd;
    end
    for i=1:nend
        if cyclicd
            im1=mod2(i-1,nd);
            ip1=mod2(i+1,nd);
        else
                im1=max(i-1,1);
                ip1=min(i+1,nd);
        end
        dn=dndt(i)*S.dt;
        dxd(i)=-dn*(yd(ip1)-yd(im1))/hypot(xd(ip1)-xd(im1),yd(ip1)-yd(im1));
        dyd(i)= dn*(xd(ip1)-xd(im1))/hypot(xd(ip1)-xd(im1),yd(ip1)-yd(im1));
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Let shoreline at shore-normal structures moves parallel to them
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %if S.struct
    %    for i=1:nend
    %        dn2(i)=dndt(i)*S.dt;
    %    end
    %    for j=1:length(str_indx)
    %        if str_indx(j)>=id1 && str_indx(j)<=id2
    %            dxd(str_indx(j))=dn2(str_indx(j)).*round(sind(thetaS(str_indx(j))),3);
    %            dyd(str_indx(j))=dn2(str_indx(j)).*cosd(thetaS(str_indx(j)));
    %        end
    %    end
    %    % dx(FS(end)+1)=dn2(FS(end)+1).*round(sind(thetaS(FS(end))),3);                                %do the same for the point after the structure
    %    % dy(FS(end)+1)=dn2(FS(end)+1).*cosd(thetaS(FS(end)));                                                %do the same for the point after the structure
    %end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Update Dune foot position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ii=1:nend
        xd(ii)=xd(ii)+dxd(ii);
        yd(ii)=yd(ii)+dyd(ii);
    end
    if cyclicd
        xd(nd+1)=xd(1);
        yd(nd+1)=yd(1);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Insert xd and yd back into x_dune,y_dune
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x_dune=[x_dune(1:id1-1),xd,x_dune(id2+1:end)];
    y_dune=[y_dune(1:id1-1),yd,y_dune(id2+1:end)];
    end
else
    x_dune=[];
    y_dune=[];
end
end