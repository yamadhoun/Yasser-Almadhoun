function [TRANSP]=get_upwindcorrection(COAST,WAVE,TRANSP,SPIT)
% function [TRANSP]=get_upwindcorrection(COAST,WAVE,TRANSP,SPIT)
%
% The upwind correction function corrects the transport to the maximum transport 
% (as computed by the angles function with an S-Phi curve made for refracted nearshore breaking waves)
% at locations where the coastline orientation is larger than the critical high angle instability angle,
% with the requirement that the previous cell was still below the critical angle (with 'ref' suffix).
%
% INPUT: 
%    COAST
%          .n           : Number of COAST.x points
%          .nq          : Number of TRANSP.QS points
%          .cyclic      : Index for COAST.cyclic coastline (0 or 1)
%    WAVE
%          .dPHItdp     : Relative angle of offshore waves with respect to the coast orientation ([Nx1] Radians)
%          .dPHIcrit    : Critical orientation of the coastline ( [Nx2] Radians)
%    TRANSP
%         option=0 : TRANSP.QS(1)=TRANSP.QSmax & Qs(2)=0; 
%         option=1 : Qs(1)=TRANSP.QSmax; Qs(2)=TRANSP.QSmax/2; Qs(3)=0; 
%         option=2 : Qs(1)=TRANSP.QSmax; Qs(2)=min(TRANSP.QS(2),TRANSP.QSmax/2); Qs(3)=0;
%          .twopoints   : Method for changes in sediment transport downdrift from the point with the upwind correction (i.e. at head of the spit) 
%          .shadowS     : Index of cells which are in the shadow zone due to coastline curvature (TRANSP.QS-points)
%          .shadowS_h   : Index of cells which are in the shadow zone due to hard structures (TRANSP.QS-points)
%          .QS          : Transport rates in grid cells [1xN] (in [m3/yr] including pores)
%          .QSmax       : Maximum transport for considered cells [1xN] (in [m3/yr] including pores)
%    SPIT
%          .spit_headwidth     (not used yet)
%
% OUTPUT:
%    TRANSP
%          .QS          : Transport rates in grid cells [1xN] (in [m3/yr] including pores)
%          .im3         : Indices of coasltine points where an upwind correction was made with positive transport (for debugging)
%          .ip3         : Indices of coasltine points where an upwind correction was made with negative transport (for debugging)
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

    debuginfo=0; % use 1, 2 or 3 for various output during debugging
    
    if isstruct(TRANSP.twopoints)
        TRANSP.twopoints=S.TRANSP.twopoints; % ensuring backward compatibility
    end

    xq=COAST.xq;
    yq=COAST.yq;
    nq=COAST.nq;
    
    TRANSP.im3=[];
    TRANSP.ip3=[];
    QS0=TRANSP.QS;
    QSmax=TRANSP.QSmax;
    if COAST.cyclic
        iirange=[1:COAST.nq];
    else
        iirange=[2:COAST.nq-1];
    end
    for cw=[1,0]
        if cw==0
            iirange=fliplr(iirange);
        end
        for i=iirange
            if COAST.cyclic
                im1=mod2(i-1,nq);
                im2=mod2(i-2,nq);
                ip1=mod2(i+1,nq);
                ip2=mod2(i+2,nq);
                idrange=mod2(i+[-5:5],nq);
            else            
                im1=max(i-1,1);
                im2=max(i-2,1);
                ip1=min(i+1,nq);
                ip2=min(i+2,nq);
                idrange=min(max(i+[-5:5],1),nq);
            end
            nq=length(TRANSP.QS);

            if length(WAVE.dPHIcrit)==1
                WAVE.dPHIcrit=repmat(WAVE.dPHIcrit,[1,length(WAVE.dPHItdp)]);
            end

            % compute relative difference between 'coast angle' and 'critical coast angle'
            % both for the current TRANSP.QS-point ('i') and the previous one (left 'im1' or right 'ip1')
            dPHIcr1(i) = mod(WAVE.dPHItdp(i)-WAVE.dPHIcrit(i)+180,360)-180;                % index showing whether currrent point is high-angle in case of positive transport direction
            dPHIcr1ref(i) = mod(WAVE.dPHItdp(im1)-WAVE.dPHIcrit(im1)+180,360)-180;         %    similar index, but then for the point just updrift of 'dPHIcr1' (for positive transport direction)
            dPHIcr2(i) = mod(WAVE.dPHItdp(i)+WAVE.dPHIcrit(i)+180,360)-180;                % index showing whether currrent point is high-angle in negative transport direction
            dPHIcr2ref(i) = mod(WAVE.dPHItdp(ip1)+WAVE.dPHIcrit(ip1)+180,360)-180;         %    similar index, but then for the point just updrift of 'dPHIcr1' (for positive transport direction)
            
            % similar variable with just the locations where dPHItdp exceeds the dPHIcrit as non zero
            % in principle this 'dPHIcr' variable may replace the 4 dPHIcr's mentioned just above, but still needs to be tested and verified first
            dPHIcr = sign(WAVE.dPHItdp).*(WAVE.dPHItdp>WAVE.dPHIcrit).*(abs(WAVE.dPHItdp)-WAVE.dPHIcrit); %dPHIcr = sign(WAVE.dPHItdp).*(WAVE.dPHItdp<-WAVE.dPHIcrit).*(abs(WAVE.dPHItdp)-WAVE.dPHIcrit);

            % compute various check variables for debugging and plotting
            dPHIrefraction=WAVE.dPHItdp-WAVE.dPHIo;
            dPHIcr0 = mod(WAVE.dPHIo-WAVE.dPHIcrit+180,360)-180;
            dQ=TRANSP.QS(2:end)-TRANSP.QS(1:end-1);
            ds=diff(COAST.s);
                         
            % plot debug info (if requested)
            if debuginfo==1
                fprintf('%3.0f, %8.6f, %8.6f (%1.0f) | %3.0f, %8.6f, %8.6f (%1.0f) | %1.0f  %1.0f  %1.0f (%1.0f)\n', ...
                i,  WAVE.dPHItdp(i),  WAVE.dPHIcrit(i),  WAVE.dPHItdp(i)>WAVE.dPHIcrit(i), ...
                im1,WAVE.dPHItdp(im1),WAVE.dPHIcrit(im2),WAVE.dPHItdp(im1)<WAVE.dPHIcrit(im2), ...
                WAVE.dPHItdp(im1)>0,~TRANSP.shadowS(im1),(isempty(TRANSP.shadowS_h)||(~isempty(TRANSP.shadowS_h)&&~TRANSP.shadowS_h(ip1)&&~TRANSP.shadowS_h(ip2))),WAVE.dPHItdp(i)>WAVE.dPHIcrit(i)&&WAVE.dPHItdp(im1)<WAVE.dPHIcrit(im2)&& WAVE.dPHItdp(im1)>0&&~TRANSP.shadowS(im1)&& (isempty(TRANSP.shadowS_h)||(~isempty(TRANSP.shadowS_h)&&~TRANSP.shadowS_h(ip1)&&~TRANSP.shadowS_h(ip2))));
            elseif debuginfo==2
                fprintf('%3.0f, %1.0f  %1.0f  %1.0f  %1.0f  %1.0f  %1.0f  %1.0f  %1.0f\n',i,WAVE.dPHItdp(i)>WAVE.dPHIcrit(i),dPHIcr1(i)>0,WAVE.dPHItdp(im1)<WAVE.dPHIcrit(im2),dPHIcr1ref(i)<=0,WAVE.dPHItdp(i)<-WAVE.dPHIcrit(i),dPHIcr2(i)<0,WAVE.dPHItdp(ip1)>-WAVE.dPHIcrit(ip2),dPHIcr2ref(i)>=0);
            end


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% UPWIND CORRECTION FOR TRANSITION POINTS TO HIGH-ANGLE WITH POSITIVE TRANSPORT DIRECTION  %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if cw==1 && dPHIcr1(i)>0 && dPHIcr1ref(i)<=0 && ... %dPHIcr2ref(i)>0 && ...
                    ~TRANSP.shadowS(im1) && COAST.clockwise==1 && ...
                    (isempty(TRANSP.shadowS_h)||(~isempty(TRANSP.shadowS_h)&&~TRANSP.shadowS_h(ip1)&&~TRANSP.shadowS_h(ip2))) 

                % correct precise position of maximum transport within vicinity (2 grid points to the left and right) of found maximum transport index i
                % id2=find(dQ(idrange(4:8))==min(dQ(idrange(4:8))),1)+4;                      % alternative method 1 : find largest change in QS (for -2/+2 grid index)
                % %id2=find(abs(dPHIcr(idrange(4:8)))==max(abs(dPHIcr(idrange(4:8)))),1)+3;   % alternative method 2 : find location with maximum difference between incoming waves and cricital angle (for -2/+2 grid index)
                % %id2=find(TRANSP.QS(idrange(4:8))==max(TRANSP.QS(idrange(4:8))),1)+3;       % alternative method 3 : find largest QS (for -2/+2 grid index)
                % i2=idrange(id2-1);
                % i1=idrange(id2);
                % ip1=idrange(id2+1);
                % ip2=idrange(id2+2);
                
                i1=i;                              % index of considered transition point (which is the first point that becomes 'high-angle')
                i2=im1;                            % index of point updrift (which is still normal 'low-angle')
                dxfactor = max(ds(i2)/ds(i1),0.5); % ratio of dx-size of grid cells (not yet used, but relevant for second order correction if grid cell sizes vary a lot)

                if TRANSP.QS(i2)<TRANSP.QS(i1) && TRANSP.QS(ip1)<TRANSP.QS(i1)         % condition where transport at i1 is already larger than for surrounding points = do nothing
                    TRANSP.QS(i1)=TRANSP.QS(i1);
                elseif TRANSP.QS(i2)>TRANSP.QS(i1) && TRANSP.QS(ip1)>TRANSP.QS(i1)     % condition where transport at i1 is smaller than for surrounding points = use mean transport
                    if TRANSP.QS(ip1)<TRANSP.QS(i2)
                        TRANSP.QS(i1)=TRANSP.QS(i2);                                   % in case of decreasing transport
                    else
                        TRANSP.QS(i1)=mean([TRANSP.QS(i2),TRANSP.QS(ip1)]);            % in case of increasing transpoirt
                    end
                elseif dPHIcr(i2)==0 && dPHIcr(ip1)>0                                  % condition where transport at i1 and just dojnwdrift (at ip1) is larger than critical angle
                    TRANSP.QS(i1)=max(QSmax([i2,i1,ip1]));           
                elseif dPHIcr(ip1)==0 && dPHIcr(i2)==0                                 % only beyond threshold for 1 cell : use max transport to stabilize coast
                    TRANSP.QS(i1)=max(TRANSP.QS([i2,i1,ip1])); 
                    if TRANSP.QS(ip1)==0
                    TRANSP.QS(i1)=max(QSmax([i2,i1,ip1]));
                    end
                else                                                                   % dPHIcr(i2)<0 && dPHIcr(ip1)>0
                    TRANSP.QS(i1)=max(QSmax([i2,i1,ip1]));                                  
                end
                %TRANSP.QS(i1)=max(min(max(QSmax(i2:i1)),TRANSP.QS(i2)),TRANSP.QS(ip1));  % alternative more complex approach accounting also for QS at updrift point (i2) and downdrift (ip1)

                % correct downdrift points in case they also exceed the dPHIcrit substantially
                maxangle=45;
                if QS0(i1)>0 && QS0(ip1)==0 && (dPHIcr(ip1)-dPHIcr(i1))>maxangle
                    TRANSP.QS(ip1)=QS0(i1);
                end
                if QS0(i1)>0 && QS0(ip1)>0 && QS0(ip2)==0 && (dPHIcr(ip2)-dPHIcr(ip1))>maxangle
                    TRANSP.QS(ip2)=QS0(ip1);
                end

                if TRANSP.twopoints==1
                    TRANSP.QS(ip1)=TRANSP.QS(ip1)+max(0,(TRANSP.QS(i1)-TRANSP.QS(ip1))/2);                %max(0,TRANSP.QS(ip1));
                elseif TRANSP.twopoints==2
                    TRANSP.QS(ip1)=TRANSP.QS(ip1)+max(0,(TRANSP.QS(i1)-TRANSP.QS(ip1))/2);                %max(0.5*TRANSP.QS(i1),TRANSP.QS(ip1));
                    TRANSP.QS(ip2)=TRANSP.QS(ip2)+max(0,(TRANSP.QS(ip1)-TRANSP.QS(ip2))/2);               %max(0,TRANSP.QS(ip2));               
                elseif TRANSP.twopoints==3                                         % NOT TESTED YET : spit width based on user defined S.spit_headwidth
                    n_spit=round((SPIT.spit_headwidth*pi)/(2*COAST.ds0));          % Calculate no. grid cells forming spit head based on half circular shape and predefined spit width
                    dist_head=linspace(1,0,n_spit);                                % Linear decay (QS_max at top, TRANSP.QS=0 at tip)
                    %dist_head=linspace(1,0,n_spit).^2;                            % Quadratic decay
                    %dist_head=exp(-linspace(0,1,n_spit));                         % Exponential decay
                    %fprintf('Twopoints 3 used, No of cells %f\n',n_spit);         % PRINT USED FOR DEBUGGING
                    for j=1:n_spit
                       if COAST.cyclic && j~=1
                           pos=mod2(i1+j-1,nq);
                       elseif j~=1
                           pos=min(i1+j-1,nq);
                       else
                           pos=i1;
                       end
                       TRANSP.QS(pos)=max(dist_head(j)*TRANSP.QS(i1),TRANSP.QS(pos));
                    end
                elseif TRANSP.twopoints==4                                         % NOT TESTED YET : spit width based on user defined S.spit_headwidth taking into account variable grid length 
                    spit_lenght=round((SPIT.spit_headwidth*pi)/(2*COAST.ds0));     % Calculate no. grid cells forming spit head based on half circular shape and predefined spit width
                    n_spit=[];
                    pos=[];
                    for j=1:spit_lenght
                       if COAST.cyclic && j~=1
                           pos(j)=mod2(i1+j-1,nq);
                           posp1=mod2(i1+j,nq);
                           n_spit(j)=hypot(x(posp1)-x(pos(j)),y(posp1)-y(pos(j)));  %determine gridcell lenght
                       elseif j~=1
                           pos(j)=min(i1+j-1,nq);
                           posp1=min(i1+j,nq);
                           n_spit(j)=hypot(x(posp1)-x(pos(j)),y(posp1)-y(pos(j)));
                       else
                           pos(j)=i1;
                           posp1=min(i1+j,nq);
                           n_spit(j)=hypot(x(posp1)-x(pos(j)),y(posp1)-y(pos(j)));
                       end
                    end
                    n_sum=cumsum(n_spit);
                    for k=1:(length(n_sum))                                                   %determine centerpoints 
                        if k==1
                            center(k)=n_sum(1)/2;                                       
                        else
                            center(k)=n_sum(k-1)+n_spit(k)/2;
                        end
                    end     
                    TRANSP.QS(pos)=max(QS_decay(center,n_spit,TRANSP.QS(i1)),TRANSP.QS(pos));                   %Apply linear decay over spit width 
                else
                    TRANSP.QS(ip1)=0;  
                end
                %elseif abs(WAVE.dPHItdp(i1))>45&&WAVE.dPHItdp(ip1)>-45&&WAVE.dPHItdp(ip1)<0
                TRANSP.ip3=[TRANSP.ip3;i1,ip1,ip2];
                
                if debuginfo==3
                    figure(99);clf;
                    hold on;
                    plot(COAST.x,COAST.y,'k-');hold on;
                    hf=fill(COAST.x,COAST.y,'y');set(hf,'FaceColor',[1 1 0.5]);
                    plot(xq([i]),yq([i]),'ks');
                    plot(xq([i2]),yq([i2]),'ks');
                    plot(xq([i1]),yq([i1]),'ko');
                    plot(xq([ip1,ip2]),yq([ip1,ip2]),'k+');
                    xlim(xq([i])+[-1000,1000]);ylim(yq([i])+[-1000,1000]);
                    title('debug plot upwind scheme')

                    fprintf('-----------------------------------------------------------------\n');
                    fprintf('name       : %8s %8s %8s %8s  \n','i2','i1','ip1','ip2');
                    fprintf('index      : %8.0f %8.0f %8.0f %8.0f  \n',[i2,i1,ip1,ip2]);
                    fprintf('QSmax      : %8.2f %8.2f %8.2f %8.2f  10^3 m^3/yr\n',QSmax([i2,i1,ip1,ip2])/10^3);
                    fprintf('QS         : %8.2f %8.2f %8.2f %8.2f  10^3 m^3/yr\n',TRANSP.QS([i2,i1,ip1,ip2])/10^3);
                    fprintf('dPHIo      : %8.3f %8.3f %8.3f %8.3f  °\n',WAVE.dPHIo([i2,i1,ip1,ip2]));
                    fprintf('dPHItdp    : %8.3f %8.3f %8.3f %8.3f  °\n',WAVE.dPHItdp([i2,i1,ip1,ip2]));
                    fprintf('dPHIo-tdp  : %8.3f %8.3f %8.3f %8.3f  °\n',WAVE.dPHIo([i2,i1,ip1,ip2])-WAVE.dPHItdp([i2,i1,ip1,ip2]));
                    fprintf('dPHIcr>THR : %8.3f %8.3f %8.3f %8.3f  °\n',dPHIcr([i2,i1,ip1,ip2]));
                    fprintf('QS(new)    : %8.2f %8.2f %8.2f %8.2f  10^3 m^3/yr\n',TRANSP.QS([i2,i1,ip1,ip2])/10^3);
                    %[QSmax(i2:ip2)/10^5;TRANSP.QS(i2:ip2)/10^5;WAVE.dPHIo(i2:ip2);WAVE.dPHItdp(i2:ip2);WAVE.dPHIo(i2:ip2)-WAVE.dPHItdp(i2:ip2);dPHIcr(i2:ip2);i2:ip2]
                end


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% UPWIND CORRECTION FOR TRANSITION POINTS TO HIGH-ANGLE WITH NEGATIVE TRANSPORT DIRECTION  %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif cw==0 && dPHIcr2(i)<0 && dPHIcr2ref(i)>=0 && ... % dPHIcr1ref(i)<0 && ...
                    ~TRANSP.shadowS(ip1) && COAST.clockwise==1 && ...
                    (isempty(TRANSP.shadowS_h)||(~isempty(TRANSP.shadowS_h)&&~TRANSP.shadowS_h(im1)&&~TRANSP.shadowS_h(im2))) 

                % correct precise position of maximum transport within vicinity (2 grid points to the left and right) of found maximum transport index i
                % id2=find(dQ(idrange(4:8))==max(dQ(idrange(4:8))),1)+3;                      % alternative method 1 : find largest change in QS (for -2/+2 grid index)
                % %id2=find(abs(dPHIcr(idrange(4:8)))==max(abs(dPHIcr(idrange(4:8)))),1)+3;   % alternative method 2 : find location with maximum difference between incoming waves and cricital angle (for -2/+2 grid index)
                % %id2=find(TRANSP.QS(idrange(4:8))==min(TRANSP.QS(idrange(4:8))),1)+3;       % alternative method 3 : find largest QS (for -2/+2 grid index)
                % i2=idrange(id2+1);
                % i1=idrange(id2);
                % im1=idrange(id2-1);
                % im2=idrange(id2-2);                
                i1=i;                              % index of considered transition point (which is the first point that becomes 'high-angle')
                i2=ip1;                            % index of point updrift (which is still normal 'low-angle')
                dxfactor = max(ds(i2)/ds(i1),0.5); % ratio of dx-size of grid cells (not yet used, but relevant for second order correction if grid cell sizes vary a lot)
                
                if TRANSP.QS(i2)>TRANSP.QS(i1) && TRANSP.QS(im1)>TRANSP.QS(i1)         % condition where transport at i1 is already larger than for surrounding points = do nothing
                    TRANSP.QS(i1)=TRANSP.QS(i1);
                elseif TRANSP.QS(i2)<TRANSP.QS(i1) && TRANSP.QS(im1)<TRANSP.QS(i1)     % condition where transport at i1 is smaller than for surrounding points = use mean transport
                    if TRANSP.QS(im1)>TRANSP.QS(i2)
                        TRANSP.QS(i1)=TRANSP.QS(i2);                                   % in case of decreasing transport
                    else
                        TRANSP.QS(i1)=mean([TRANSP.QS(i2),TRANSP.QS(im1)]);            % in case of increasing transpoirt
                    end
                elseif dPHIcr(i2)==0 && dPHIcr(im1)<0                                  % condition where transport at i1 and just dojnwdrift (at im1) is larger than critical angle
                    TRANSP.QS(i1)=-1*max(QSmax([i2,i1,im1]));                                 
                elseif dPHIcr(i2)==0 && dPHIcr(im1)==0                                 % only beyond threshold for 1 cell : use max transport to stabilize coast
                    TRANSP.QS(i1)=min(TRANSP.QS([i2,i1,im1])); 
                    if TRANSP.QS(im1)==0
                    TRANSP.QS(i1)=-1*max(QSmax([i2,i1,im1]));
                    end
                else                                                                   % dPHIcr(i2)>0 && dPHIcr(im1)<0
                    TRANSP.QS(i1)=-1*max(QSmax([i2,i1,im1]));  
                    %TRANSP.QS(i1)=min(max(-max(QSmax(i1:i2)),TRANSP.QS(i2)),TRANSP.QS(im1));  % alternative more complex approach accounting also for QS at updrift point (i2) and downdrift (ip1)
                end  
                % TRANSP.QS(i1)=TRANSP.QS(i1)+(1-max(TRANSP.QS(im1)/-QSmax(im1),0))*(-QSmax(i1)-TRANSP.QS(i1)); % gradual approach using fraction of difference between QSmax and QS(i)

                % correct downdrift points in case they also exceed the dPHIcrit substantially
                maxangle=45;
                if QS0(i1)<0 && QS0(im1)==0 && -(dPHIcr(im1)-dPHIcr(i1))>maxangle
                    TRANSP.QS(im1)=QS0(i1);
                end
                if QS0(i1)<0 && QS0(im1)<0 && QS0(im2)==0 && -(dPHIcr(im2)-dPHIcr(im1))>maxangle
                    TRANSP.QS(im2)=QS0(im1);
                end

                if TRANSP.twopoints==1 
                    TRANSP.QS(im1)=TRANSP.QS(im1)+min(0,(TRANSP.QS(i1)-TRANSP.QS(im1))/2);                % min(0,TRANSP.QS(im1));
                elseif TRANSP.twopoints==2
                    TRANSP.QS(im1)=TRANSP.QS(im1)+min(0,(TRANSP.QS(i1)-TRANSP.QS(im1))/2);                %min(0.5*TRANSP.QS(i1),TRANSP.QS(im1));
                    TRANSP.QS(im2)=TRANSP.QS(im2)+min(0,(TRANSP.QS(im1)-TRANSP.QS(im2))/2);               %min(0,TRANSP.QS(im2));
                elseif TRANSP.twopoints==3                                         % NOT TESTED YET : spit width based on user defined S.spit_headwidth
                    n_spit=round((SPIT.spit_headwidth*pi)/(2*COAST.ds0));          % Calculate no. grid cells forming spit head based on half circular shape and predefined spit width
                    dist_head=linspace(1,0,n_spit);                                % Linear decay (QS_max at top, TRANSP.QS=0 at tip)
                    for j=1:n_spit
                       if COAST.cyclic && j~=1
                           pos=mod2(i1-j+1,nq);
                       elseif j~=1
                           pos=min(i1-j+1,nq);
                       else
                           pos=i1;
                       end
                       TRANSP.QS(pos)=-1*max(dist_head(j)*TRANSP.QS(i1),TRANSP.QS(pos));
                    end
                elseif TRANSP.twopoints==4                                         % NOT TESTED YET : spit width based on user defined S.spit_headwidth taking into account variable grid length 
                    spit_lenght=round((SPIT.spit_headwidth*pi)/(2*COAST.ds0));     % Calculate no. grid cells forming spit head based on half circular shape and predefined spit width
                    n_spit=[];
                    pos=[];
                    for j=1:spit_lenght
                       if COAST.cyclic && j~=1
                           pos(j)=mod2(i1-j+1,nq);
                           posp1=mod2(i1-j,nq);
                           n_spit(j)=hypot(x(posp1)-x(pos(j)),y(posp1)-y(pos(j)));      %determine gridcell lenght
                       elseif j~=1
                           pos(j)=min(i1-j+1,nq);
                           posp1=min(i1-j,nq);
                           n_spit(j)=hypot(x(posp1)-x(pos(j)),y(posp1)-y(pos(j)));
                       else
                           pos(j)=i1;
                           posp1=min(i1-j,nq);
                           n_spit(j)=hypot(x(posp1)-x(pos(j)),y(posp1)-y(pos(j)));
                       end
                    end
                    n_sum=cumsum(n_spit);
                    for k=1:(length(n_sum))                                                   %determine centerpoints 
                        if k==1
                            center(k)=n_sum(1)/2;                                       
                        else
                            center(k)=n_sum(k-1)+n_spit(k)/2;
                        end
                    end     
                    TRANSP.QS(pos)=-1*max(QS_decay(center,n_spit,TRANSP.QS(i1)),TRANSP.QS(pos));                   %Apply linear decay over spit width 
                else
                    TRANSP.QS(im1)=0;  
                end
                TRANSP.im3=[TRANSP.im3;i1,im1,im2];

                if debuginfo==3
                    figure(99);clf;
                    hold on;
                    plot(COAST.x,COAST.y,'k-');hold on;
                    hf=fill(COAST.x,COAST.y,'y');set(hf,'FaceColor',[1 1 0.5]);
                    hold on;hp3=plot(xq([im1,im2,i1]),yq([im1,im2,i1]),'r+');
                    hp1=plot(xq([i2]),yq([i2]),'rs');           
                    hp2=plot(xq([i1]),yq([i1]),'ro');
                    xlim(xq([i])+[-1000,1000]);ylim(yq([i])+[-1000,1000]);
                    title('debug plot upwind scheme')

                    fprintf('-----------------------------------------------------------------\n');
                    fprintf('name       : %8s %8s %8s %8s  \n','im2','im1','i1','i2');
                    fprintf('index      : %8.0f %8.0f %8.0f %8.0f  \n',[im1,im2,i1,i2]);
                    fprintf('QSmax      : %8.2f %8.2f %8.2f %8.2f  10^3 m^3/yr\n',-QSmax([im1,im2,i1,i2])/10^3);
                    fprintf('QS         : %8.2f %8.2f %8.2f %8.2f  10^3 m^3/yr\n',QS0([im1,im2,i1,i2])/10^3);
                    fprintf('dPHIo      : %8.3f %8.3f %8.3f %8.3f  °\n',WAVE.dPHIo([im1,im2,i1,i2]));
                    fprintf('dPHItdp    : %8.3f %8.3f %8.3f %8.3f  °\n',WAVE.dPHItdp([im1,im2,i1,i2]));
                    fprintf('dPHIo-tdp  : %8.3f %8.3f %8.3f %8.3f  °\n',WAVE.dPHIo([im1,im2,i1,i2])-WAVE.dPHItdp([im1,im2,i1,i2]));
                    fprintf('dPHIcr>THR : %8.3f %8.3f %8.3f %8.3f  °\n',dPHIcr([im1,im2,i1,i2]));
                    fprintf('QS(new)    : %8.2f %8.2f %8.2f %8.2f  10^3 m^3/yr\n',TRANSP.QS([im1,im2,i1,i2])/10^3);
                    %[-QSmax(im2:i2)/10^5;TRANSP.QS(im2:i2)/10^5;WAVE.dPHIo(im2:i2);WAVE.dPHItdp(im2:i2);WAVE.dPHIo(im2:i2)-WAVE.dPHItdp(im2:i2);dPHIcr(im2:i2);im2:i2]
                end
            end
        end
    end
    
    % debug value
    TRANSP.debug.QS2=TRANSP.QS;
end