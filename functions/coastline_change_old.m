function [COAST,GROYNE]=coastline_change(COAST,TRANSP,STRUC,GROYNE,TIME,NOUR,CC)
% function [COAST.x,COAST.y,S,WAVE.indxw,COAST.dx,COAST.dSds,sgroyne,xgroyne,ygroyne]=coastline_change(COAST.ds0,COAST.h0,NOUR.growth,TRANSP.boundary_cond_start,TRANSP.boundary_cond_end,FORMAT.SLplot,COAST.s,TRANSP.QS,COAST.x,COAST.y,STRUC,GROYNE,GROYNE.QS,FORMAT.tplot,TIME.tnow,TIME.tend,TIME.adt,NOUR.rate_m3_per_m_yr,TIME.it,TIME.dt,WAVE.WVC,WAVE.indxw,COAST.n_mc,i_mc,BATHY)
% 
% updates the coastline at each time step based on transport gradients and nourishment volumes
% modifications to describe coastline behaviour next to groynes
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

    eps=1e-1;
    idgroyne=[];

    %% spatial steps
    n=COAST.n;
    nq=COAST.nq;
    nend=n;

    %% Spatial discretisation
    x0=COAST.x;y0=COAST.y;
    ndev=zeros(1,n);
    COAST.dSds=zeros(1,n);
    COAST.dndt=zeros(1,n);
    COAST.dx=zeros(1,n);
    COAST.dy=zeros(1,n);
    COAST.ds=zeros(1,n);
    
    %% Determine dominant bypass transport
    for ig=1:GROYNE.n
        [QSm,imx]=max(abs(GROYNE.QS(ig,:)));
        GROYNE.QS(ig,:)=GROYNE.QS(ig,imx);
    end
    
    %% Coastline treatment if groyne at beginning of coastline section
    %  Note: need to do this before x,y are updated
    x1 = nan;
    y1 = nan;
    
    if ~isempty(GROYNE.x)
        %[idgroyne]=find(hypot(GROYNE.x(:,2)-COAST.x(1),GROYNE.y(:,2)-COAST.y(1))<eps);  % find connected groyne
        [idgroyne]=find(GROYNE.idcoast(:,2)==COAST.i_mc);
    end

    if ~isempty(idgroyne)
        if GROYNE.QS(idgroyne,2)>0
            %% bypass over shadow zone
            %ishadow=find(abs(TRANSP.QS)>eps,1,'first')-1;
            %ishadow=find(TRANSP.QS>eps,1,'first')-1;
            ishadow=find(TRANSP.QS>eps,1,'first');
            if ishadow>0
                TRANSP.QS(1:ishadow)=GROYNE.QS(idgroyne,2);
            end
        end
        ist=GROYNE.strucnum(idgroyne);
        [~    ,xhard,~   ,~ , ~ ] = get_one_polygon( STRUC.s_hard,STRUC.x_hard,ist );
        [shard,yhard,~   ,~ , ~ ] = get_one_polygon( STRUC.s_hard,STRUC.y_hard,ist );
        COAST.ds(1)=hypot(COAST.x(2)-COAST.x(1),COAST.y(2)-COAST.y(1))/2.0;
        COAST.dSds(1)=(TRANSP.QS(1)-GROYNE.QS(idgroyne,2))/COAST.ds(1);
        COAST.dndt(1)=-COAST.dSds(1)/COAST.h0;
        dn=(COAST.dndt(1)+(NOUR.growth*NOUR.rate_m3_per_m_yr(1)/COAST.h0)-CC.SLRo/COAST.tanbeta)*TIME.dt;
        dsg=-dn;
        shard_tip=shard(GROYNE.tipindx(idgroyne,2));
        GROYNE.s(idgroyne,2)=max(GROYNE.s(idgroyne,2)+dsg,shard_tip);
        GROYNE.x(idgroyne,2)=interp1(shard,xhard,GROYNE.s(idgroyne,2));
        GROYNE.y(idgroyne,2)=interp1(shard,yhard,GROYNE.s(idgroyne,2));
        if ~isnan(GROYNE.x(idgroyne,2))
            x1=GROYNE.x(idgroyne,2);
            y1=GROYNE.y(idgroyne,2);
        end
    end

    %% Coastline treatment if groyne at end of coastline section
    xend = nan;
    yend = nan;
    if ~isempty(GROYNE.x)
        %[idgroyne]=find(hypot(GROYNE.x(:,1)-COAST.x(end),GROYNE.y(:,1)-COAST.y(end))<eps); % find connected groyne
        [idgroyne]=find(GROYNE.idcoast(:,1)==COAST.i_mc);
    end
    if ~isempty(idgroyne)
        if GROYNE.QS(idgroyne,1)<0
            %% bypass over shadow zone
            %ishadow=find(abs(TRANSP.QS)>eps,1,'last')+1;
            %ishadow=find(TRANSP.QS<-eps,1,'last')+1;
            ishadow=find(TRANSP.QS<-eps,1,'last');
            if ishadow<length(TRANSP.QS)
                TRANSP.QS(ishadow:end)=GROYNE.QS(idgroyne,1);
            end
        end
        ist=GROYNE.strucnum(idgroyne);
        [~    ,xhard,~   ,~ , ~ ] = get_one_polygon( STRUC.s_hard,STRUC.x_hard,ist );
        [shard,yhard,~   ,~ , ~ ] = get_one_polygon( STRUC.s_hard,STRUC.y_hard,ist );
        COAST.ds(end)=hypot(COAST.x(end)-COAST.x(end-1),COAST.y(end)-COAST.y(end-1))/2;
        COAST.dSds(end)=(GROYNE.QS(idgroyne,1)-TRANSP.QS(end))/COAST.ds(end);
        COAST.dndt(end)=-COAST.dSds(end)/COAST.h0;                          % dndt in [m/yr] = alongshore gradient in sediment transport in [m3/yr/m] divided by 'active height' [m]
        dn=(COAST.dndt(end)+(NOUR.growth*NOUR.rate_m3_per_m_yr(end)/COAST.h0)-CC.SLRo/COAST.tanbeta)*TIME.dt;      % dndt in [m/yr] + nourishment in [m3/m/yr] divided by the active height [m]
        dsg=dn;
        shard_max=shard(GROYNE.tipindx(idgroyne,1));
        if shard_max==0
            warning('shard_max=0, coastline_change')
        end
        %GROYNE.s(idgroyne,1)=min(GROYNE.s(idgroyne,1)+dsg,shard_max);
        GROYNE.s(idgroyne,1)=GROYNE.s(idgroyne,1)+dsg;
        if GROYNE.s(idgroyne,1)>shard_max
            disp('GROYNE.s(idgroyne,1)>shard_max')
            GROYNE.s(idgroyne,1)=shard_max;
        end
        GROYNE.x(idgroyne,1)=interp1(shard,xhard,GROYNE.s(idgroyne,1));
        GROYNE.y(idgroyne,1)=interp1(shard,yhard,GROYNE.s(idgroyne,1));
        if ~isnan(GROYNE.x(idgroyne,1))
            xend=GROYNE.x(idgroyne,1);
            yend=GROYNE.y(idgroyne,1);
        end
    end
   
    %% COMPUTE RATE OF COASTLINE CHANGE
    % dSds : rate of volumetric change in [m3/yr per meter coastline length]
    % dndt : rate of coastline change in [m/yr]
    for i=1:n
        if COAST.cyclic
            im1=mod2(i-1,n);
            ip1=mod2(i+1,n);
            im1q=mod2(i-1,nq);
            ip1q=mod2(i,nq);
        else
            im1=max(i-1,1);
            ip1=min(i+1,n);
            im1q=max(i-1,1);
            ip1q=min(i,nq);
        end
        
        COAST.ds(i)=hypot(COAST.x(ip1)-COAST.x(im1),COAST.y(ip1)-COAST.y(im1));
        COAST.dSds(i)=(TRANSP.QS(ip1q)-TRANSP.QS(im1q))*2/COAST.ds(i); %,COAST.ds0);  
        COAST.dndt(i)=-COAST.dSds(i)/COAST.h0;                                      
    end
    
    %% ADD BOUNDARY CONDITIONS
    if strcmpi(TRANSP.boundary_condition_start,'Periodic') && ~COAST.cyclic
        % periodic b.c.
        COAST.dndt(1)=COAST.dndt(end-1);
    end
    if strcmpi(TRANSP.boundary_condition_end,'Periodic') && ~COAST.cyclic
        % periodic b.c.
        COAST.dndt(end)=COAST.dndt(2);
    end
    if strcmpi(TRANSP.boundary_condition_start,'Closed') && ~COAST.cyclic
        % zero transport b.c.
        COAST.dSds(1)=TRANSP.QS(1)/COAST.ds(1);
        COAST.dndt(1)=-COAST.dSds(1)/COAST.h0;
    end
    if strcmpi(TRANSP.boundary_condition_end,'Closed') && ~COAST.cyclic
        % zero transport b.c.
        COAST.dSds(nend)=-TRANSP.QS(nend-1)/COAST.ds(nend-1);
        COAST.dndt(nend)=-COAST.dSds(nend)/COAST.h0;
    end
    if strcmpi(TRANSP.boundary_condition_start,'Neumann') && ~COAST.cyclic
        % zero transport b.c.
        COAST.dSds(1)=0;
        COAST.dndt(1)=-COAST.dSds(1)/COAST.h0;
    end
    if strcmpi(TRANSP.boundary_condition_end,'Neumann') && ~COAST.cyclic
        % zero transport b.c.
        COAST.dSds(nend)=0;
        COAST.dndt(nend)=-COAST.dSds(nend)/COAST.h0;
    end

    for i=1:n
        if COAST.cyclic
            im1=mod2(i-1,n);
            ip1=mod2(i+1,n);
        else
            im1=max(i-1,1);
            ip1=min(i+1,n);
        end
        dn=(COAST.dndt(i)+(NOUR.rate_m3_per_m_yr(i)/COAST.h0)-CC.SLRo/COAST.tanbeta)*TIME.dt;
        COAST.dx(i)=-dn*(COAST.y(ip1)-COAST.y(im1))/COAST.ds(i); %,COAST.ds0);
        COAST.dy(i)= dn*(COAST.x(ip1)-COAST.x(im1))/COAST.ds(i); %,COAST.ds0);
    end

    %% Update coastline positions, all points
    for i=1:n
        COAST.x(i)=x0(i)+COAST.dx(i);
        COAST.y(i)=y0(i)+COAST.dy(i);
    end
    if COAST.cyclic
        COAST.x(n)=COAST.x(1);
        COAST.y(n)=COAST.y(1);
    end

    %% Overwrite coastline position if groyne at begin of coastline section
    if ~isnan(x1)
        COAST.x(1)=x1;
        COAST.y(1)=y1;
    end
    %% Overwrite coastline position if groyne at end of coastline section
    if ~isnan(xend)
        COAST.x(n)=xend;
        COAST.y(n)=yend;
    end

%     figure(100);clf;
%     plot(x0,y0,'b.-');hold on;plot(COAST.x,COAST.y,'r-');
%     fprintf('%1.1f ',COAST.dndt);
%     fprintf('\n');
%     a=1;
    % A2a = area(polyshape([x0(:),y0(:)]));
    % A2b = area(polyshape([COAST.x(:),COAST.y(:)]));fprintf('A2a=%12.0f, A2b=%12.0f => diff=%12.0f \n',A2a,A2b,A2b-A2a);
end