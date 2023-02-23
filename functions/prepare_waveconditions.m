function [WAVE]=prepare_waveconditions(S,TIME)
% function [WAVE]=prepare_waveconditions(S,TIME)
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

    fprintf('  Prepare wave conditions \n');
    
    indxw=[];
    WC=struct;
    WVC=struct;
    
    WAVE=struct;
    WAVE.WC=struct;
    WAVE.WVC=struct;
    WAVE.ddeep=S.ddeep;
    WAVE.dnearshore=S.dnearshore;
    WAVE.gamma=S.gamma;
    WAVE.wave_interaction=S.wave_interaction;
    WAVE.tide_interaction=S.tide_interaction;
    WAVE.wavefile=S.wavefile;
    WAVE.surf_width=S.surf_width;
    WAVE.surf_width_w=S.surf_width_w;
    WAVE.WVCfile=S.WVCfile;
    WAVE.spread=S.spread;                                                             % wave spreading [?] (wave_dir from range:  S.phiw0 +/- 0.5*S.spread)
    WAVE.Waveclimfile=S.Waveclimfile;                                                 % wave climate file
    WAVE.Wavematfile=S.Wavematfile;
    WAVE.Wavecorr=S.Wavecorr;   
    WAVE.interpolationmethod=S.interpolationmethod;  
    WAVE.diffraction=S.diffraction;    
    WAVE.spacevaryingwave=0;    
    try WAVE.DA=S.DA;end
    WAVE.ccSLR=S.ccSLR;               % climate impacted rise in sea level (SLR) [Nx2] with 'time in datenum format' and 'sea level with respect to initial situation' (rates per year are computed automatically) % S.tanbeta is used as 'slope angle'.
    WAVE.ccDIR=S.ccDIR;               % climate impacted change in wave direction [Nx2] with 'time in datenum format' and 'relative change in wave direction w.r.t. initial situation' as a # degrees (rates per year are computed automatically) 
    WAVE.ccHS=S.ccHS;                 % climate impacted change in wave height [Nx2] with 'time in datenum format' and 'relative change in wave height w.r.t. initial situation' as a non-dimensionless multiplicationfactor (rates per year are computed automatically) 
    
    if ~isempty(WAVE.WVCfile)
        if ~isfield(WAVE,'WVCtimestep')
            WAVE.WVCtimestep = TIME.dt*365;
        end  
        if ischar(WAVE.WVCfile)
            WAVE.WVCfile = {WAVE.WVCfile};
        end       
        if (size(WAVE.WVCfile,2)>1 && size(WAVE.WVCfile,1)==1)
            WAVE.WVCfile=WAVE.WVCfile';
        end
        if isfield(WAVE,'WVCtime')
            WAVE.WVCtime=S.WVCtime; 
        end

        %% CONVERT FILES THAT REFERENCE MULTIPLE WAVE CLIMATE FILES/LOCATIONS
        % convert file format if a WVT/WVC file is used as input with reference to a large number of WVT/WVC files 
        if (strcmpi(WAVE.WVCfile{1}(end-3:end),'.WVT') ...
           || strcmpi(WAVE.WVCfile{1}(end-3:end),'.WVC')) && length(WAVE.WVCfile)==1
            try 
                WVCraw=load(WAVE.WVCfile{1});
            catch 
                WVCfid=fopen(WAVE.WVCfile{1},'r');
                WAVE.WVCfile = {};
                gg=0;
                while ~feof(WVCfid)
                    WVCline={};
                    try
                        gg=gg+1;
                        WVCline=fgetl(WVCfid);
                    end
                    WVCinfo = textscan(WVCline,'%s %f %f');
                    [dirnm,filnm,extnm]=fileparts(WVCinfo{1}{1});
                    if isempty(dirnm)
                       dirnm='.'; 
                    end    
                    WAVE.WVCfile{gg,1} = [dirnm,filesep,filnm,extnm];
                    WAVE.WVCfile{gg,2} = WVCinfo{2};
                    WAVE.WVCfile{gg,3} = WVCinfo{3};
                end
            end
        end
        
        % convert file format if a GKL-file is used as input with reference to a large number of RAY files 
        if strcmpi(WAVE.WVCfile{1}(end-3:end),'.GKL') && length(WAVE.WVCfile)==1
            GKLdata=readGKL(WAVE.WVCfile{1});
            [dirnm,filnm,extnm]=fileparts(WAVE.WVCfile{1});
            if isempty(dirnm)
                dirnm='.';
            end
            WAVE.WVCfile = {};
            for gg=1:length(GKLdata.x)
                WAVE.WVCfile{gg,1} = [dirnm,filesep,GKLdata.ray_file{gg},'.ray'];
                WAVE.WVCfile{gg,2} = GKLdata.x(gg);
                WAVE.WVCfile{gg,3} = GKLdata.y(gg);
            end
        end
             
        %% LOOP OVER WAVE TIME SERIES POINTS IN SPACE
        kk1=1;
        warning off
        WVC(size(WAVE.WVCfile,1))=struct;
        WAVE.spacevaryingwave=size(WAVE.WVCfile,1)>1;
        for kk=1:size(WAVE.WVCfile,1)
            fprintf('   reading file : %s\n',WAVE.WVCfile{kk});
            
            % READ NC-FILES FROM SnapWave
            if strcmpi(WAVE.WVCfile{kk}(end-2:end),'.nc')
                lat=ncread(WAVE.WVCfile{kk},'latitude');
                long=ncread(WAVE.WVCfile{kk},'longitude');

                Raduis = 6378137.0;
                x = Raduis * cosd(lat) .* cosd(long);
                y = Raduis * cosd(lat) .* sind(long);

                hm0=ncread(WAVE.WVCfile{kk},'depth');
                Hs=ncread(WAVE.WVCfile{kk},'Hs');
                tp=ncread(WAVE.WVCfile{kk},'Tp');
                wd=ncread(WAVE.WVCfile{kk},'Dm');

                timeunits=ncreadatt(WAVE.WVCfile{kk},'time','units');
                reftime=datenum(timeunits(15:end)); %reftime=datenum(1976,01,01);
                time=reftime+ncread(WAVE.WVCfile{kk},'time')/24/60/60; 
                for kk2=1:length(x)
                    WVC(kk1).x=x(kk2);
                    WVC(kk1).y=y(kk2);
                    WVC(kk1).timenum=time;
                    if length(x)==size(hm0,1)
                        WVC(kk1).Hs=squeeze(hm0(kk2,:));
                        WVC(kk1).Tp=squeeze(tp(kk2,:));
                        WVC(kk1).Dir=squeeze(wd(kk2,:));
                    else
                        WVC(kk1).Hs=squeeze(hm0(:,kk2));
                        WVC(kk1).Tp=squeeze(tp(:,kk2));
                        WVC(kk1).Dir=squeeze(wd(:,kk2));
                    end
                    kk1=kk1+1;
                end
                WAVE.spacevaryingwave=length(WVC)>1;
            
            % READ WAVE TIME-SERIES DATA
            elseif strcmpi(WAVE.WVCfile{kk}(end-3:end),'.WVT') 
                WVCraw=load(WAVE.WVCfile{kk});warning off;
                if length(num2str(WVCraw(1,1)))==12
                    WVC(kk).timenum=datenum([num2str(WVCraw(:,1))],'yyyymmddHHMM');
                elseif length(num2str(WVCraw(1,1)))<12
                    WVC(kk).timenum=WVCraw(:,1);
                    %WVC(kk).timenum=datenum([num2str(WVCraw(:,1))],'yyyymmdd');
                else
                    WVC(kk).timenum=WVCraw(:,1);
                    %WVC(kk).timenum=datenum([num2str(WVCraw(:,1))],'yyyymmddHHMMSS');       
                 end
                WVC(kk).Hs=interpNANs(WVCraw(:,2));
                WVC(kk).Tp=interpNANs(WVCraw(:,3));
                WVC(kk).Dir=interpNANsDIR(WVCraw(:,4));
            
            % READ WAVE TIME-SERIES DATA
            elseif strcmpi(WAVE.WVCfile{kk}(end-3:end),'.asc')
                WVCraw=load(WAVE.WVCfile{kk});warning off;
                WVC(kk).timenum=datenum([num2str(WVCraw(:,1),'%08.0f'),num2str(WVCraw(:,2),'%06.0f')],'yyyymmddHHMMSS');
                WVC(kk).Hs=interpNANs(WVCraw(:,3));
                WVC(kk).Tp=interpNANs(WVCraw(:,4));
                WVC(kk).Dir=interpNANsDIR(WVCraw(:,5));
                
            % WAVE CLIMATE
            elseif strcmpi(WAVE.WVCfile{kk}(end-3:end),'.WVC')
                WVCraw=load(WAVE.WVCfile{kk});warning off;
                WVC(kk).timenum=[];
                WVC(kk).Hs=interpNANs(WVCraw(:,1))';
                WVC(kk).Tp=interpNANs(WVCraw(:,2))';
                WVC(kk).Dir=interpNANsDIR(WVCraw(:,3))';
                if size(WVCraw,2)>3
                   WVC(kk).Prob=interpNANs(WVCraw(:,4))';    
                   
                   % scale to 1.0
                   WVC(kk).Prob=WVC(kk).Prob/sum(WVC(kk).Prob);
                end
                
            % TIME-SERIES OR STATIC RAY files
            elseif strcmpi(WAVE.WVCfile{kk}(end-3:end),'.RAY')
                RAYdata=readRAY(WAVE.WVCfile{kk,1});warning off;
                WVC(kk).time     = [];  
                WVC(kk).timenum  = TIME.timenum0;    % use as fall-back option the starttime of the simulation as starttime of the WVC-file
                if isfield(WAVE,'WVCtime') 
                    WVC(kk).timestr  = WAVE.WVCtime;
                    WVC(kk).timenum  = datenum(WAVE.WVCtime); % use time specified per file
                end
                try  
                    WVC(kk).time = RAYdata.time;     % import/read time (only in case of time-series ray-file)
                end
                if ~isempty(RAYdata.time)
                    WVC(kk).timenum  = WVC(kk).timenum+RAYdata.time/24;
                end
                WVC(kk).PHIequi  = RAYdata.Cequi;                  % angle of equilibrium also accounting for the net currents which are in QSoffset (=computeEQUI(WVC(kk).equi,WVC(kk).c1,WVC(kk).c2,WVC(kk).hoek,WVC(kk).QSoffset) 
                WVC(kk).c1       = RAYdata.c1*10^6;     
                WVC(kk).c2       = RAYdata.c2; 
                WVC(kk).h0       = RAYdata.h0;    
                WVC(kk).fshape   = RAYdata.fshape; 
                WVC(kk).PHIf     = RAYdata.hoek;
                WVC(kk).QSoffset = zeros(size(WVC(kk).PHIequi));
                if ~isempty(RAYdata.QSoffset)
                    WVC(kk).QSoffset = RAYdata.QSoffset; % timeseries ray with tide offset
                end
                WVC(kk).Hs        = (abs(WVC(kk).c1).^0.5)/200.0;    % use a proxy     Hs=sqrt(c1)/200
                WVC(kk).Tp        = WVC(kk).Hs*2+3;                  % use a proxy     Tp=Hs*2+3
                WVC(kk).Dir       = RAYdata.hoek-RAYdata.equi;       % angle of equilibrium accounting only for the wave part (i.e. hoek-equi)              % use a proxy     Dir=Cequi
            
            % TIME-SERIES FROM WAVE-MAT-FILE FOR DATA ASSIMILATION OR COUPLING
            elseif  strcmpi(WAVE.WVCfile{kk}(end-3:end),'.mat')
                WVCraw=load(WAVE.WVCfile{kk});warning off;
                if WAVE.DA==1
                    WVC(kk).timenum=WVCraw.WaveDA_record(:,1);
                    WVC(kk).Hs=interpNANs(WVCraw.WaveDA_record(:,2));
                    WVC(kk).Tp=interpNANs(WVCraw.WaveDA_record(:,3));
                    WVC(kk).Dir=interpNANsDIR(WVCraw.WaveDA_record(:,4));
                else
                    WVC(kk).timenum=datenum([num2str(WVCraw.WVCn(:,1))],'yyyymmdd');
                    WVC(kk).Hs=interpNANs(WVCraw.WVCn(:,2));
                    WVC(kk).Tp=interpNANs(WVCraw.WVCn(:,3));
                    WVC(kk).Dir=interpNANsDIR(WVCraw.WVCn(:,4));
                end
            
            % IF NO SUITBALE INPUT IS PROVIDED
            else
                fprintf('Warning : Wave condition file cannot be read. Please check the formatting! \n')
            end
            
            % ADD LOCATION X AND Y
            if ~strcmpi(WAVE.WVCfile{kk}(end-2:end),'.nc')
                WVC(kk).x = kk;
                WVC(kk).y = kk;
                try
                    WVC(kk).x = WAVE.WVCfile{kk,2};
                    WVC(kk).y = WAVE.WVCfile{kk,3};
                catch
                    if kk>1
                        fprintf('%s\n',' prepare_wave_conditions :: Warning : please specify an X and Y location of each of the WVCfiles (i.e. second and third column of S.WVCfile) \n');
                    end
                end
            end
        end
        
        %% AGGREGATE TIME-SERIES DATA
        % method of speeding up the simulation by computing average conditions for time step periods
        if isfield(WVC(1),'timenum') && ~isempty(WVC(1).timenum)
            for kk=1:length(WVC)
                DT=diff(WVC(kk).timenum);
                DT0=median(DT);
                factor=max(ceil(WAVE.WVCtimestep/DT0),1);
                DT1=factor*DT0;
                ttend=min(ceil(length(WVC(kk).timenum)/factor)+1,length(WVC(kk).timenum));
                
                if factor>2
                    fldvc = {'Dir' 'PHIequi' 'PHIf'};
                    fldsc = {'Hs' 'Tp' 'c1' 'c2' 'h0' 'fshape' 'QSoffset'};
                    sfac = [2.0 1.0 1.0 1.0 1.0 1.0 1.0]; %scaling factors for scalars
                    
                    WVCtmp=struct;
                    for tt=1:ttend
                        ttrange=unique(min(max((tt-1)*factor+[-floor((factor-1)/2):ceil((factor-1)/2)],1),length(WVC(kk).timenum)));
                        WVCtmp.timenum(tt,1) = WVC(kk).timenum(1)+DT1*(tt-1);
                        
                        % aggregate scalars over period 'ttrange'
                        for jj=1:length(fldsc)
                            if isfield(WVC(kk),fldsc{jj})
                                WVCtmp.(fldsc{jj})(tt,1) = (mean(WVC(kk).(fldsc{jj})(ttrange).^sfac(jj))).^(1./sfac(jj));
                            end
                        end
                        % aggregate vectors over period 'ttrange'
                        for jj=1:length(fldvc)
                            if isfield(WVC(kk),fldvc{jj})
                                HSfac = sfac(1);
                                sdir0=sind(WVC(kk).(fldvc{jj})(ttrange)); sdir0 = sdir0(:);
                                cdir0=cosd(WVC(kk).(fldvc{jj})(ttrange)); cdir0 = cdir0(:);
                                sdir1=mean(sdir0.*WVC(kk).Hs(ttrange).^HSfac ./ sum(WVC(kk).Hs(ttrange).^HSfac) );
                                cdir1=mean(cdir0.*WVC(kk).Hs(ttrange).^HSfac ./ sum(WVC(kk).Hs(ttrange).^HSfac) );
                                WVCtmp.(fldvc{jj})(tt,1) = mod(atan2d(sdir1,cdir1),360);
                            end
                        end
                    end
                    
                    WVC(kk).timenum      = WVCtmp.timenum;
                    for jj=1:length(fldsc)
                        if isfield(WVC(kk),fldsc{jj})
                            WVC(kk).(fldsc{jj}) = WVCtmp.(fldsc{jj});
                        end
                    end
                    for jj=1:length(fldvc)
                        if isfield(WVC(kk),fldvc{jj})
                            WVC(kk).(fldvc{jj}) = WVCtmp.(fldvc{jj});
                        end
                    end
                end
                
                % additional index of timepoints indw introduced by Ahmed (=time with respect to start of simulation in days)
                WVC(kk).indw=[];
                for iw=1:length(WVC(kk).timenum)
                    WVC(kk).indw(iw)=TIME.timenum0-WVC(kk).timenum(iw);
                end
                WVC(kk).indxw=0;
                if ~isempty(find(WVC(kk).indw<=0))
                    WVC(kk).indxw=find(WVC(kk).indw==max(WVC(kk).indw(WVC(kk).indw<=0)));
                end
                
                % check if timeseries covers valid time range (after model start / 'timenum0')
                if ~isempty(WVC(kk).timenum)
                    if WVC(kk).timenum(1) > TIME.timenum0 %Warning Message
                        warning(['prepare_wave_conditions :: Error: The model start time is not covered by the wave data for timeseries file ',num2str(kk)]);
                    end
                end
            end
        end

    %% SPECIAL CASE : WAVEclimfile used as a proxy for a spatially uniform climate 
    elseif ~isempty(WAVE.Waveclimfile)
        WCraw=load(WAVE.Waveclimfile);
        WC=struct;
        WC.Hs=WCraw(:,1);
        WC.Tp=WCraw(:,2);
        WC.dir=WCraw(:,3);  %+WAVE.Wavecorr;
        if size(WCraw,2)==4
            WC.Prob= WCraw(:,4)/sum(WCraw(:,4));  
        end
    end

    %% export wave variables
    WAVE.WVC=WVC;
    WAVE.WC=WC;
end
