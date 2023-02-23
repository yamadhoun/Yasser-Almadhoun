function [qs,qw,qss,qww]=dune_evolution(S,x_dune,y_dune,x_dune0,y_dune0,x_mc,y_mc,x_mc0,y_mc0,i_mc,Dfelev,phiwnd,SWL,cyclic,Hso,str_indx,shadowD,shadowD_h,philoc)
% function [qs,qw,qss,qww]=dune_evolution(S,x_dune,y_dune,x_dune0,y_dune0,x_mc,y_mc,x_mc0,y_mc0,i_mc,Dfelev,phiwnd,SWL,cyclic,Hso,str_indx,shadowD,shadowD_h,philoc)
%
% computes the cumulative distance along a line
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ESTIMATE DUNE FOOT CHANGE DUE TO INTERACTION 
%% BETWEEN WAVE ACTION AND AEOLIAN TRANSPORT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if S.dune
    %empirical parameters
    kf=0.02;                                % Friction coefficient
    Cs=5e-7;                                % Impact coefficient
    d50r=2.5e-4;                            % Median reference grain size
    rhoa=1.225;                             % Air density
    Aw=0.1;                                 % Coefficient
    Kw=4.2;                                 % Imperical coefficient
    k=0.41;                                 % Von Karman's coefficient
    segmaw=0.1;                             % empirical factor
    phid=zeros(size(x_dune));
    [ xc,yc,~,~,~ ] = get_one_polygon( x_mc,y_mc,i_mc );               % In case of multiple coastline sections--> use each section seperately
    [ xc0,yc0,~,~,~ ] = get_one_polygon( x_mc0,y_mc0,i_mc );      % In case of multiple coastline sections--> use each section seperately
    qs=zeros(size(xc));
    qw=zeros(size(xc));
    len=hypot(max(xc)-min(xc),max(yc)-min(yc));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find number of sections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nansd=find(isnan(x_dune));
    n_dune=length(nansd)+1;
    for i_dune=1:n_dune
        [ xd,yd,~,id1,id2] = get_one_polygon( x_dune,y_dune,i_dune ); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estimate dune foot angle w.r.t  x-axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nd0=length(xd)-1;
        phid0=zeros(size(xd));
        nend=nd0+1;
        for i=1:nend
            im1=max(i-1,1);
            ip1=min(i+1,nd0+1);
            phid0(i)=360-atan2d(yd(ip1)-yd(im1),xd(ip1)-xd(im1));
        end
        phid(id1:id2)=phid0;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estimate local wave and wind angles of each dune foot grid point == No. dune foot grid points 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for im=1:length(xc)-1
        phic(im)=360-atan2d(yc(im+1)-yc(im),xc(im+1)-xc(im));
    end
    phicwloc=atan2d(sind(phic-phiwnd),cosd(phic-phiwnd));
    sc(1)=0;                 
    for ix=2:length(xc)
        sc(ix)=sc(ix-1)+hypot(xc(ix)-xc(ix-1),yc(ix)-yc(ix-1));               % accum. distances of coastline grid points
    end
    xdune=x_dune;
    ydune=y_dune;
    xdune(nansd)=[];
    ydune(nansd)=[];
    for ji=1:length(xdune)
        Xss=[xdune(ji)-len*sind(phid(ji)),xdune(ji)+len*sind(phid(ji))]; 
        Yss=[ydune(ji)-len*cosd(phid(ji)),ydune(ji)+len*cosd(phid(ji))];
        [xs1,ys1]=intersectfn(xc,yc,Xss,Yss);  
        if ~isempty(xs1)
            [~, indx_sc] = min(hypot(xc-xs1(1),yc-ys1(1)));
            ds2=hypot(xc(indx_sc+1)-xs1(1),yc(indx_sc+1)-ys1(1));
            if ds2<=sc(indx_sc+1)-sc(indx_sc)
                phidloc(ji)=philoc(indx_sc);
                phidwloc(ji)=phicwloc(indx_sc);
            else
                phidloc(ji)=philoc(indx_sc-1);
                phidwloc(ji)=phicwloc(indx_sc-1);
            end
        else
            phidloc(ji)=0;
            phidwloc(ji)=0;
        end
    end
    if n_dune>1
        phidloc=[phidloc(1:nansd-1),0,phidloc(nansd:end)];           % needed to be updated for more than one nan
        phidwloc=[phidwloc(1:nansd-1),0,phidwloc(nansd:end)];        % needed to be updated for more than one nan
        phidwloc(abs(phidwloc)>90)=90;                                   % set wind angle to 90 for 90 <angles< -90 
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estimate width of the berm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Wberm0=zeros(size(x_dune));
        Wberm=zeros(size(x_dune));
        for iii=1:length(x_dune)
             Xh=[x_dune(iii)-len*sind(phid(iii)),x_dune(iii)+len*sind(phid(iii))]; 
             Yh=[y_dune(iii)-len*cosd(phid(iii)),y_dune(iii)+len*cosd(phid(iii))];
             [xx1,yy1]=intersectfn(xc,yc,Xh,Yh);
                if ~isempty(xx1)
                    for ij=1:length(xx1)
                        all_dis(ij) = hypot(x_dune(iii)-xx1(ij),y_dune(iii)-yy1(ij));
                    end
                    Wberm(iii)=min(all_dis);                         % existing  berm width     
                    all_dis=[];
                end
        end
        for lll=1:length(x_dune0)
            Xs0=[x_dune0(lll)-len*sind(phid(lll)),x_dune0(lll)+len*sind(phid(lll))]; 
            Ys0=[y_dune0(lll)-len*cosd(phid(lll)),y_dune0(lll)+len*cosd(phid(lll))];
            [xx2,yy2]=intersectfn([xc0;yc0],[Xs0;Ys0]);
            if ~isempty(xx2)
                Wberm0(lll)=hypot((x_dune0(lll)-xx2(1)),(y_dune0(lll)-yy2(1)));                   % initial berm width        
            end
        end
        if Wberm(1)==0;  Wberm(1)=Wberm(2); end
        if Wberm(end)==0; Wberm(end)=Wberm(end-1);end
        if Wberm0(1)==0;Wberm0(1)=Wberm0(2);end
        if Wberm0(end)==0;Wberm0(end)=Wberm0(end-1);end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estimate beach slope
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for ii=1:length(x_dune)
              slope(ii)=Dfelev(ii)/Wberm(ii);
              if slope(ii)>0.5                         % when wberm bocomes too small dramatically when stime step is too long
                  slope(ii)=0.5;                        % to avoid having very high erosion rates due to very steep slope
              end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% transport rate of eroded sediment from the dune to the beach qs [m2/yr] per m width of dune foot
%% Estimate Adjusted Run-up height & deepwater wave length & Mean wave period
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Hso1 =Hso.*sqrt(cosd(phidloc));                          % wave height indicates onshore energy flux [m]
        Hso1(abs(phidloc)>90)=0;                       % Ignore offshore directed waves 
        L0=1.56*S.tper.^2;                                     %Deepwater wave length [m]
        R=slope.*sqrt(Hso1.*L0);                           % Run-up Height [m]
        TH=mean(S.tper);
        R1=zeros(size(x_dune));
        for jj=1:length(R)
            if R(jj)+SWL>Dfelev(jj)
                R1(jj)=R(jj)*exp(-2*kf*Wberm(jj))+Dfelev(jj)*(1-exp(-2*kf*Wberm(jj)));    % Adjusted run-up height due to berm friction
                qss(jj)=4*Cs*(((R1(jj)-Dfelev(jj)+SWL).^2)/TH)*60*60*24*365;                    % Transport rate due to wave attack [m2/yr]/m'
            else
                qss(jj)=0;
            end
        end
        qss(nansd)=0;
        qss(shadowD)=0;                                           % Set dune transport to zero due to dune-dune shadowing
        qss(shadowD_h)=0;                                       % Set dune transport to zero due to being sheltered behind structures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Aeolian transport rate from beach to dune
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        uc=Aw*sqrt((S.rhos-rhoa)*S.g*S.d50/rhoa);                                     %wind critical shear velocity (m/s)
        uc=uc.*ones(1,length(x_dune));
        z0=S.d50/30;
        us=k*S.uz/log(S.z/z0);                                                                         %wind shear velocity (m/s)  
        us=us.*ones(1,length(x_dune));
        for ii=1:length(us)
            if us(ii)>uc(ii)
                mwe(ii)=Kw*sqrt(S.d50/d50r)*rhoa*us(ii)^2*(us(ii)-uc(ii))/S.g;                       %potential aeolian transport rate (Kg/m/s)
            else
                mwe(ii)=0;
            end
        end
        % Equlibruim transport rate [m2/yr] 
        qwe=60*60*24*365*mwe*(1-S.porosity)/S.rhos;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estimation of the dry beach width
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nd=length(x_dune)-1;
        for i=1:nend
            if cyclic
                im1=mod2(i-1,nd);
                ip1=mod2(i+1,nd);
            else
                im1=max(i-1,1);
                ip1=min(i+1,nd+1);
            end
        end
        for kk=1:length(x_dune)
            if (R(kk)+SWL)<Dfelev(kk)
                Bdry(kk)=((1-((R(kk)+SWL)./Dfelev(kk))).*Wberm(kk));
            else
                Bdry(kk)=0;
            end
        end
        FL=Bdry./cosd(phidwloc);
        qw0=qwe.*(1- (exp(-segmaw*FL)));
        qww=qw0.*cosd(phidwloc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Boundary Conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        qss(1)=0;qww(1)=0;
        qss(end)=0;qww(end)=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Interaction between Coastline and dunefoot 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for jjj=2:length(x_dune)-1
            Xs=[x_dune(jjj)-len*sind(phid(jjj)),x_dune(jjj)+len*sind(phid(jjj))]; 
            Ys=[y_dune(jjj)-len*cosd(phid(jjj)),y_dune(jjj)+len*cosd(phid(jjj))];
            [xx1,yy1]=intersectfn(xc,yc,Xs,Ys);
            if ~isempty(xx1)
                [~, indx] = min(hypot(xc-xx1,yc-yy1));
            end
            qs(indx)=qss(jjj);
            qw(indx)=qww(jjj);
        end
        if S.struct
            qs(str_indx)=0;qw(str_indx)=0;   % to keep the downdrift coastline grid point in place as it is always sheltered by the structure
        end
else
    qs=[];
    qw=[];
    qss=[];
    qww=[];
end
end