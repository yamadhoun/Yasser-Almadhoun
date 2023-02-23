function [STRUC,WAVE] = wave_diffraction(STRUC,COAST,WAVE);
%% wave diffraction following Roelvink approach for angles, Kamphuis for Kd
%
% INPUT:
%    STRUC
%        .x_hard      : x-coordinate of hard structures
%        .y_hard      : y-coordinate of hard structures
%    COAST
%        .x           : x-coordinate of coastline (only current section)
%        .y           : y-coordinate of coastline (only current section)
%        .n           : Number of coastline section points
%    WAVE
%        .PHItdp      : Wave direction at the toe of the dynamic profile
%        .HStdp       : Wave height at the toe of the dynamic profile [m]
%
% OUTPUT:
%     STRUC
%         .xtip,ytip  : coordinates of diffraction points
%         .Htip       : Wave height at structure tips
%         .dirtip     : Wave directions at structure tips
%     WAVE
%        .PHItdp      : Wave direction after diffraction
%        .HStdp       : Wave height after diffraction
%
%% Copyright notice
%   --------------------------------------------------------------------
%   Copyright (C) 2020 IHE Delft & Deltares
%
%       Dano Roelvink, Ahmed Elghandour
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

    rotfac = 0.8;
    omegat = -20;
    kdform = 'Roelvink';
    %% Get input
    x_hard = STRUC.x_hard;
    y_hard = STRUC.y_hard;
    x      = COAST.x;
    if length(COAST.x)<3
       disp('This is too short a coastline segment, good sir.') 
       pause 
    end
    y      = COAST.y;
    xS     = 0.5*(x(1:end-1)+x(2:end));
    yS     = 0.5*(y(1:end-1)+y(2:end));
    nS     = length(xS);
    PHItdp = WAVE.PHItdp;
    HStdp  = WAVE.HStdp;
    
    %% Determine diffraction points xtip,ytip for all structures
    [ x_struc,y_struc,n_hard,~,~ ] = get_one_polygon( x_hard,y_hard,1);
    for j=1:n_hard
        [ x_struc,y_struc,n_hard,~,~ ] = get_one_polygon( x_hard,y_hard,j );
        [ wetstr,wetstr,n_hard,~,~ ] = get_one_polygon( STRUC.wetstr_mc,STRUC.wetstr_mc,j );
        % Project structure coordinates on line in wave direction
        [~,closest]=min(hypot(mean(x_struc)-xS,mean(y_struc)-yS));
        PHItip=PHItdp(closest);
        proj=(x_struc-x_struc(1))*cosd(PHItip)-(y_struc-y_struc(1))*sind(PHItip);
        % The minimum and maximum distances to this line are the diffraction
        % points
        proj(wetstr==0)=1e10;
        [~,itip1]=min(proj);
        proj(wetstr==0)=-1e10;
        [~,itip2]=max(proj);
        xtip(1,j)=x_struc(itip1);
        ytip(1,j)=y_struc(itip1);
        xtip(2,j)=x_struc(itip2);
        ytip(2,j)=y_struc(itip2);
        for tip=1:2
            [~,closest]=min(hypot(xtip(tip,j)-xS,ytip(tip,j)-yS));
            Htip(tip,j)=HStdp(closest);
            dirtip(tip,j)=PHItdp(closest);
        end
    end
    
    %% Compute angle to diffraction point for each transport point,
    %% relative to incident wave angle
    % phitip1 and phitip2 are angles from coast transport point to diffraction
    % points [n_hard by nS]
    phitip1=atan2d(xtip(1,:)'-xS,ytip(1,:)'-yS); % matrix [n_hard by nS]
    phitip2=atan2d(xtip(2,:)'-xS,ytip(2,:)'-yS);
    % omega1 and omega2 angles within diffraction zone, positive inwards
    omega1=repmat(dirtip(1,:)',1,nS)-phitip1;     % matrix [n_hard by nS]
    omega2=phitip2-repmat(dirtip(2,:)',1,nS);
    omega1=mod(omega1+180,360)-180;
    omega2=mod(omega2+180,360)-180;
    
    %% Determine which points are in which influence zone
    in=omega1>omegat&omega2>omegat;              % matrix [n_hard by nS]
    
    %% om is the change of direction due to diffraction
    om=zeros(size(HStdp));
    % om1, om2 are wave directions relative to incident direction due to each
    % diffraction point, rotated by % rotfac (=2 according to Hurst et al,
    % but 0.8 according to comparison with % wave-resolving XBeach model.
    % Both om1 and om2 are zero at edge of influence zone.
    om1=min(max(0,(omega1-omegat)*rotfac),90);
    om2=-min(max(0,(omega2-omegat)*rotfac),90);
    
    %% Diffraction coefficients according to simple formula
    % kdform= 'Roelvink' or 'Kamphuis'
    kd1=wave_diffraction_coeff(omega1,kdform);
    kd2=wave_diffraction_coeff(omega2,kdform);
    
    %% sumin is the 1d array denoting the influence areas of structures:
    % 0 = no structure, 1 = in lee of one structure, 2 = between two structures
    % and influenced by both.
    sumin=sum(in,1);
    for j=1:n_hard
        inn=in(j,:); % 1d mask for influence of structure j
        
        %% Check if diffraction points can be 'seen' from coast transport points
        % If not, set kd to 0. epsx and epsy are to avoid counting
        % intersections in offshore parts of structure with finite width
        for i=find(inn)
            epsx=0.2*(xtip(1,j)-xS(i));
            epsy=0.2*(ytip(1,j)-yS(i));
            [xx_hard1,yy_hard1]=intersectfn([xS(i) xtip(1,j)-epsx],[yS(i) ytip(1,j)-epsy],x_hard,y_hard);
            epsxc=0.001*(xtip(1,j)-xS(i));
            epsyc=0.001*(ytip(1,j)-yS(i));
            if length(find(inn))>1
                [xx_coast1,yy_coast1]=intersectfn([xS(i)+epsxc xtip(1,j)],[yS(i)+epsyc ytip(1,j)],x(find(inn)),y(find(inn)));
                if ~isempty(xx_hard1) | ~isempty(xx_coast1)
                    kd1(j,i)=0;
                end
            end
            epsx=0.2*(xtip(2,j)-xS(i));
            epsy=0.2*(ytip(2,j)-yS(i));
            [xx_hard2,yy_hard2]=intersectfn([xS(i) xtip(2,j)-epsx],[yS(i) ytip(2,j)-epsy],x_hard,y_hard);
            epsxc=0.001*(xtip(2,j)-xS(i));
            epsyc=0.001*(ytip(2,j)-yS(i));
            if length(find(inn))>1
                [xx_coast2,yy_coast2]=intersectfn([xS(i)+epsxc xtip(2,j)],[yS(i)+epsyc ytip(2,j)],x(find(inn)),y(find(inn)));
                if ~isempty(xx_hard2) | ~isempty(xx_coast2)
                    kd2(j,i)=0;
                end
            end
        end
        %% Combine wave height and direction from both sides
        HStdp(inn)=sqrt((kd1(j,inn)*Htip(1,j)).^2+(kd2(j,inn)*Htip(2,j)).^2);
        om(inn)=atan2d((kd1(j,inn)*Htip(1,j)).^2.*sind(om1(j,inn)) ...
                      +(kd2(j,inn)*Htip(2,j)).^2.*sind(om2(j,inn)), ...
                       (kd1(j,inn)*Htip(1,j)).^2.*cosd(om1(j,inn)) ...
                      +(kd2(j,inn)*Htip(2,j)).^2.*cosd(om2(j,inn)));
    end
    %% Points in the overlapping influence areas of two structures
    % multiply kd, do weighted averaging of directions
    for i=find(sumin==2)
        kd=1;
        sinkd2=0;
        coskd2=0;
        for j=find(in(:,i))'
            [kdj,ind]=max([kd1(j,i),kd2(j,i)]);
            kd=kd*kdj;
            Ht=0;
            if ind==1
                sinkd2=sinkd2+kdj^2*sind(om1(j,i));
                coskd2=coskd2+kdj^2*cosd(om1(j,i));
                Ht=Ht+Htip(1,j);
            else
                sinkd2=sinkd2+kdj^2*sind(om2(j,i));
                coskd2=coskd2+kdj^2*cosd(om2(j,i));
                Ht=Ht+Htip(2,j);
            end
        end
        HStdp(i)=kd*Ht/2;
        om(i)=atan2d(sinkd2,coskd2);
    end
    
    %% Output
    STRUC.xtip  = xtip;
    STRUC.ytip  = ytip;
    STRUC.Htip  = Htip;
    STRUC.xp    = xS;    
    STRUC.yp    = yS;
    WAVE.PHItdp = PHItdp-om;  % minus because om is cartesian, PHItdp nautical
    WAVE.HStdp  = HStdp;
    WAVE.dirtip = dirtip;
    
end
    