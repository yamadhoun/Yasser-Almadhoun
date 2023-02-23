function [V,FORMAT,TIME]=plot_coast(CHANNEL,STRUC,COAST,CLIFF,WAVE,TIME,TRANSP,FORMAT,V)
% function [V,FORMAT,TIME]=plot_coast(CHANNEL,STRUC,COAST,WAVE,TIME,FORMAT,V)
%
% INPUT:
%   STRUC
%      .x_hard
%      .y_hard
%   COAST
%      .x_mc0
%      .y_mc0
%      .PHIf
%      .x_mc
%      .y_mc
%   WAVE
%      .HSo_mc
%      .PHIo_mc
%      .HSbr_mc
%      .PHIbr_mc
%      .WVC
%   TIME
%      .it
%      .tnow
%      .dt
%      .tnext
%   TRANSP
%      .QS_mc
%   FORMAT
%      .tplot
%      .tsl
%      .plot_time
%      .xmax
%      .ymax
%      .xmin
%      .ymin
%      .nmax
%      .xp_mc
%      .yp_mc
%      .xlimits
%      .ylimits
%      .XYoffset
%      .SLplot
%      .fignryear
%      .outputdir
%      .video
%      .plotinterval;
%      .usefill;
%      .LDBplot;
%      .llocation
%      .ld
%   V
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
%   This library is free software: you can redistribute TIME.it and/or
%   modify TIME.it under the terms of the GNU Lesser General Public
%   License as published by the Free Software Foundation, either
%   version 2.1 of the License, or (at your option) any later version.
%
%   This library is distributed in the hope that TIME.it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with this library. If not, see <http://www.gnu.org/licenses
%   --------------------------------------------------------------------

    qvrscale=1.5;
    eps=-1d-10;
%     if isempty(WAVE.HStdp_mc)
%         WAVE.HStdp_mc=WAVE.HSo_mc.*0;
%         WAVE.PHIbr_mc=WAVE.PHIo_mc;
%     end
%     if isempty(WAVE.HSbr_mc)
%         WAVE.HSbr_mc=WAVE.HSo_mc.*0;
%         WAVE.PHIbr_mc=WAVE.PHIo_mc;
%     end
%     WAVE.HSbr_mc(isnan(WAVE.HSbr_mc))=0;
%     IDnotnan=~isnan(WAVE.HSo_mc);
%     HS_mean = median(WAVE.HSo_mc(WAVE.HSo_mc-mean(WAVE.HSo_mc(IDnotnan))>eps));       % necessary for spatially uniform waves
%     PHI_mean = median(WAVE.PHIo_mc(WAVE.HSo_mc-mean(WAVE.HSo_mc(IDnotnan))>eps));
%     HStdp_mean = median(WAVE.HStdp_mc(WAVE.HSo_mc-mean(WAVE.HSo_mc(IDnotnan))>eps));
%     PHItdp_mean = median(WAVE.PHItdp_mc(WAVE.HSo_mc-mean(WAVE.HSo_mc(IDnotnan))>eps));
%     HSbr_mean = median(WAVE.HSbr_mc(WAVE.HSo_mc-mean(WAVE.HSo_mc(IDnotnan))>eps));
%     PHIbr_mean = median(WAVE.PHIbr_mc(WAVE.HSo_mc-mean(WAVE.HSo_mc(IDnotnan))>eps));
    
    if ~exist(fullfile(pwd,FORMAT.outputdir),'dir')
        mkdir(fullfile(pwd,FORMAT.outputdir));
    end
    
    if mod(TIME.it,FORMAT.plotinterval)==0
        %iplot=iplot+1;
        if 1
            %% FILL THE COASTLINE
            [FORMAT.xmax,FORMAT.ymax,FORMAT.xmin,FORMAT.ymin,FORMAT.nmax ] = plot_fill_sections(COAST.x_mc,COAST.y_mc,FORMAT.ld,TIME.it,FORMAT.xmax,FORMAT.ymax,FORMAT.xmin,FORMAT.ymin,FORMAT.nmax,COAST.ds0,FORMAT.usefill);
            hold on;
            plot(STRUC.x_hard,STRUC.y_hard,'k','linewidth',2);
            plot(COAST.x_mc0,COAST.y_mc0,'b','linewidth',1); %2
            plot(CLIFF.x,CLIFF.y,'r','linewidth',2); %2
            if ~isempty(CHANNEL.xr_mc') 
                plot(CHANNEL.xr_mc,CHANNEL.yr_mc,'b.-')
                plot(CHANNEL.x_inlet,CHANNEL.y_inlet,'*')
            end
            
            %% PLOT TRANSPORTS (optional)
            if FORMAT.plotQS==1
                try
                    hsc1=scatter(COAST.x_mc(1:end-1),COAST.y_mc(1:end-1),4,TRANSP.QS_mc);
                    
                    set(hsc1,'markerFaceColor','flat');
                    set(gca,'Clim',[-800000,800000]);hold on;
                end
            end

            %% NOT YET USED PLOTS FOR DUNES -> requires WIND and DUNE data
            %if S.dune && S.qplot
            % [qwave,qwind,innt,time,step]=plot_qrates_wave_wind(S,x_dune,y_dune,x_trans,y_trans,n_trans,qss,qww,innt,TIME.tnow,qwave,qwind,time,step);   % <- new function by AE&MG
            %end
            %if S.dune && S.bermw_plot
            % [BWplot, bermW,int,BWplot2]=plot_berm_width(S,COAST.x_mc,COAST.y_mc,x_dune,y_dune,xl,yl,n_trans,BWplot,BWplot2,bermW,TIME.tnow,int,x_trans,y_trans);    % <- new function by AE&MG
            %end
            %if S.CLplot
            % [CLplot,int,CLplot2,ds_cl]=plot_coastline_transect(S,COAST.x_mc,COAST.y_mc,COAST.x_mc0,COAST.y_mc0,n_trans,CLplot,CLplot2,TIME.tnow,iint,x_trans,y_trans,ds_cl);    % <- new function by AE&MG
            %end
            

            %% PLOT OFFSHORE WAVE HEIGHT TEXT
%             if isempty(FORMAT.XYwave)
%                 FORMAT.XYwave = [FORMAT.XYoffset(1)+FORMAT.xlimits(1)+(sind(PHI_mean)+1)/2*diff(FORMAT.xlimits),FORMAT.XYoffset(2)+FORMAT.ylimits(1)+(cosd(PHI_mean)+1)/2*diff(FORMAT.ylimits),diff(FORMAT.xlimits)/12];
%             end

            %% PLOT WAVE QUIVER & OFFSHORE WAVE HEIGHT TEXT
%             PHI2=[];
%             PHItdp2=[];
%             PHIbr2=[];
%             PHIf2=[];
%             HS2=[];
%             HStdp2=[];
%             HSbr2=[];
%             xx=[];
%             yy=[];  
%             if isempty(WAVE.WVC) || WAVE.spacevaryingwave==0
%                 arx=[FORMAT.XYwave(1)-FORMAT.XYoffset(1)];
%                 ary=[FORMAT.XYwave(2)-FORMAT.XYoffset(2)];
%                 PHI2=PHI_mean;
%                 PHItdp2=PHItdp_mean;
%                 PHIbr2=PHIbr_mean;
%                 HS2=HS_mean;
%                 HStdp2=HStdp_mean;
%                 HSbr2=HSbr_mean;
%                 if ~isempty(COAST.PHIf) 
%                     PHIf2=median(COAST.PHIf);
%                 end
%                 xx=min(max(FORMAT.XYwave(1)-FORMAT.XYoffset(1)-diff(FORMAT.xlimits)/24,min(FORMAT.xlimits)),max(FORMAT.xlimits));
%                 yy=min(max(FORMAT.XYwave(2)-FORMAT.XYoffset(2)+diff(FORMAT.ylimits)/24,min(FORMAT.ylimits)),max(FORMAT.ylimits));
%             else %if ~isempty(COAST.x_mc)
%                 WVCid=zeros(length(WAVE.WVC),1);
%                 arx=zeros(length(WAVE.WVC),1);
%                 ary=zeros(length(WAVE.WVC),1);
%                 PHI2=zeros(length(WAVE.WVC),1);
%                 PHItdp2=zeros(length(WAVE.WVC),1);
%                 PHIbr2=zeros(length(WAVE.WVC),1);
%                 HS2=zeros(length(WAVE.WVC),1);
%                 HStdp2=zeros(length(WAVE.WVC),1);
%                 HSbr2=zeros(length(WAVE.WVC),1);
%                 xx=zeros(length(WAVE.WVC),1);
%                 yy=zeros(length(WAVE.WVC),1);
%                 PHIf2=zeros(length(WAVE.WVC),1);
%                 for gg=1:length(WAVE.WVC)
%                     try
%                         dist=((COAST.x_mc-WAVE.WVC(gg).x).^2+(COAST.y_mc-WAVE.WVC(gg).y).^2).^0.5;
%                         WVCid(gg)=min(find(dist==min(dist),1),length(WAVE.PHIo_mc));
%                         arx(gg)=[WAVE.WVC(gg).x-FORMAT.XYoffset(1)];
%                         ary(gg)=[WAVE.WVC(gg).y-FORMAT.XYoffset(2)];
%                         PHI2(gg)=WAVE.PHIo_mc(WVCid(gg));
%                         PHItdp2(gg)=WAVE.PHItdp_mc(WVCid(gg));
%                         PHIbr2(gg)=WAVE.PHIbr_mc(WVCid(gg));
%                         HS2(gg)=WAVE.HSo_mc(WVCid(gg));
%                         HStdp2(gg)=WAVE.HStdp_mc(WVCid(gg));
%                         HSbr2(gg)=WAVE.HSbr_mc(WVCid(gg));
%                         xx(gg)=arx(gg);
%                         yy(gg)=ary(gg);
%                         if ~isempty(COAST.PHIf)
%                             PHIf2(gg)=COAST.PHIf(min(WVCid(gg),length(COAST.PHIf)));
%                             xx(gg)=arx(gg) + diff(FORMAT.xlimits)/35 .* sind(PHIf2(gg));
%                             yy(gg)=ary(gg) + diff(FORMAT.ylimits)/35 .* cosd(PHIf2(gg));
%                         end
%                     end
%                 end
%             end
%             for gg=1:length(PHI2)
%                 factor = HS2(gg)/max(HS2);
%                 hqvr1=quiver(arx(gg),ary(gg),FORMAT.XYwave(3).*cosd(3/2*180-PHI2(gg))*factor,FORMAT.XYwave(3).*sind(3/2*180-PHI2(gg))*factor,qvrscale);hold on;
%                 set(hqvr1,'linewidth',2,'Color','k','AutoScale','off','AutoScaleFactor',1.2,'MaxHeadSize',1);
%                 hqvr2=quiver(arx(gg),ary(gg),FORMAT.XYwave(3).*cosd(3/2*180-PHItdp2(gg))*HStdp2(gg)/HS2(gg)*factor,FORMAT.XYwave(3).*sind(3/2*180-PHItdp2(gg))*HStdp2(gg)/HS2(gg)*factor,qvrscale);hold on;
%                 set(hqvr2,'linewidth',2,'Color',[0 0 0.8],'AutoScale','off','AutoScaleFactor',1.2,'MaxHeadSize',1);
%                 hqvr3=quiver(arx(gg),ary(gg),FORMAT.XYwave(3).*cosd(3/2*180-PHIbr2(gg))*HSbr2(gg)/HS2(gg)*factor,FORMAT.XYwave(3).*sind(3/2*180-PHIbr2(gg))*HSbr2(gg)/HS2(gg)*factor,qvrscale);hold on;
%                 set(hqvr3,'linewidth',2,'Color',[0.8 0 0],'AutoScale','off','AutoScaleFactor',1.2,'MaxHeadSize',1);
%                 htxt=text(xx(gg),yy(gg),['',num2str(HSbr2(gg),'%2.1f'),'m']);
%                 set(htxt,'FontSize',9,'HorizontalAlignment','Center','VerticalAlignment','Middle','Color',[0.8 0 0]);
%             end

            %% TRY PLOT LOWER SHOREFACE ORIENTATION
%             if ~isempty(PHIf2) && ~strcmpi(TRANSP.trform,'KAMP') && ~strcmpi(TRANSP.trform,'CERC') && ~strcmpi(TRANSP.trform,'CERC2')
%                 for gg=1:length(arx)
%                     dx=FORMAT.XYwave(3).*cosd(PHIf2(gg))/3.*[1,-1];
%                     dy=FORMAT.XYwave(3).*sind(PHIf2(gg))/3.*[-1,1];
%                     hl=plot(arx(gg)+dx,ary(gg)+dy);hold on;
%                     set(hl,'linewidth',1,'Color',[0 0.5 0]);
%                 end
%             end

            %% PLOT REFERENCE LINES AND LEGEND
            if ~isempty(FORMAT.LDBplot)
                for mm=1:size(FORMAT.LDBplot,1)
                    try
                        LDBplotval = get_landboundary(FORMAT.LDBplot{mm,1});
                    catch
                        LDBplotval = load(FORMAT.LDBplot{mm,1});
                    end
%                     hp6(mm)=plot(LDBplotval(:,1)-FORMAT.XYoffset(1),LDBplotval(:,2)-FORMAT.XYoffset(2),FORMAT.LDBplot{mm,3},'linewidth',1.5);
                end
%                 hleg = legend(hp6,FORMAT.LDBplot(:,2)','Location',FORMAT.llocation);
%                 set(hleg,'Box','off','Color','None');
            end

            %% PLOT SHORELINES AT SPECIFIC DATES
            if ~isempty(FORMAT.SLplot)&& TIME.dt==(FORMAT.tplot-TIME.tnow)/365
                try
                FORMAT.xp_mc{FORMAT.tsl,:}=COAST.x_mc;
                FORMAT.yp_mc{FORMAT.tsl,:}=COAST.y_mc;
                FORMAT.tsl=FORMAT.tsl+1;
                FORMAT.tplot=FORMAT.plot_time(FORMAT.tsl);
                end
            else%if isempty(FORMAT.SLplot)
                FORMAT.xp_mc=FORMAT.xp_mc;
                FORMAT.yp_mc=FORMAT.yp_mc;
            end
            if ~isempty(FORMAT.SLplot) && FORMAT.tsl>1
                if ~isempty(FORMAT.LDBplot)
                    LName=FORMAT.LDBplot(:,2);
                else
                    LName={};
                    mm=0;
                end

                for sl=1:size(FORMAT.xp_mc,1)
%                     hp6(sl+mm)=plot(FORMAT.xp_mc{sl,:},FORMAT.yp_mc{sl,:},FORMAT.SLplot{sl,3},'linewidth',1.5);
                    LName(end+1)=FORMAT.SLplot(sl,2);
                end
%                 hleg = legend(hp6,LName','Location',FORMAT.llocation);
%                 set(hleg,'Box','off','Color','None');
            elseif isempty(FORMAT.SLplot) && ~isempty(FORMAT.LDBplot)
%                 hleg = legend(hp6,FORMAT.LDBplot(:,2)','Location',FORMAT.llocation);
%                 set(hleg,'Box','off','Color','None');
            end
        end
        
        %% FORMATTING THE PLOT
        hold off;
        ax = get(FORMAT.mainfighandle,'CurrentAxes');
        set(ax,'DataAspectRatio',[1 1 1],'DataAspectRatioMode','manual','PlotBoxAspectRatio',[3 4 4],'PlotBoxAspectRatioMode','manual'); % manual axis equal
        xlim(FORMAT.xlimits);
        ylim(FORMAT.ylimits);
        set(ax,'XtickLabel',num2str(get(ax,'Xtick')'/1000,'%2.1f'));
        set(ax,'YtickLabel',num2str(get(ax,'Ytick')'/1000,'%2.1f'));
        xlabel(ax,'Easting [km]','HandleVisibility','off');
        ylabel(ax,'Northing [km]','HandleVisibility','off');
        date=datevec(TIME.tnow);
        title(ax,num2str(date(1:3)));
        if ~FORMAT.fastplot
           drawnow;    % 10% of plottime
        end
        
        %% video
        if FORMAT.video==1 && STRUC.diffraction==0
            V(TIME.it+1)=getframe(FORMAT.mainfighandle);
        end
        
        %% images
        if TIME.tnow>=TIME.tnext
            fname=[num2str(round((TIME.it+1)),'%04.0f')]; 
            if ~FORMAT.fastplot    
               print(FORMAT.mainfighandle,fullfile(pwd,FORMAT.outputdir,[fname '.jpg']),'-djpeg','-r300');
            else
               F=getframe(gcf);
               imwrite(F.cdata,fullfile(pwd,FORMAT.outputdir,[fname '.png']),'png');
            end
            vii=length(V)+1;
            V(vii)=getframe(FORMAT.mainfighandle);
            TIME.tnext=TIME.tnext+365./FORMAT.fignryear;
        end
        FORMAT.xyt(TIME.it+1).COAST.x_mc=COAST.x_mc;
        FORMAT.xyt(TIME.it+1).COAST.y_mc=COAST.y_mc;
        
    end
end
