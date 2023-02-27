close all; clc;
%% ===============================subplots=================================
figure('WindowState','maximized'); 
subplot(2,2,1)
plot(O.timenum(1,:),O.dndt_cliff(1,:),'linewidth',2);
datetick
set(gca,'FontSize',12)
title('Cliff Retreat Over Time in Years 2022-2052 Using TWL Model')
xlabel('Time [year]')
ylabel('Cliff Retreat [m/yr]')
%%
subplot(2,2,2)
plot(O.BeachWidth(1,:),O.dndt_cliff(1,:),'linewidth',2);
set(gca,'FontSize',12)
title('Cliff Retreat Versus Beach Width in Years 2022-2052 Using TWL Model')
xlabel('Beach Width [m]')
ylabel('Cliff Retreat [m/yr]')
%%
subplot(2,2,3)
plot(O.timenum(1,:),O.BeachWidth(1,:),'linewidth',2);
datetick
set(gca,'FontSize',12)
title('Beach Width Over Time in Years 2022-2052 Using TWL Model')
xlabel('Time [year]')
ylabel('Beach Width [m]')
%%
subplot(2,2,4)
plot(O.timenum(1,:),O.cliffposition(1,:),'linewidth',2); 
datetick
set(gca,'FontSize',12)
title('Cliff Position Over Time in Years 2022-2052 Using TWL Model')
xlabel('Time [year]')
ylabel('Cliff Position [m]')


%% ============================cliff/coast position========================
figure('WindowState','maximized'); hold on
plot(O.timenum(1,:),O.cliffposition(1,:),'linewidth',2); 
plot(O.timenum(1,:),O.coastposition(1,:),'linewidth',2); 
datetick
set(gca,'FontSize',18)
title('Cliff Position Over Time in Years 2022-2052 Using TWL Model')
xlabel('Time [year]')
ylabel('Cliff Position [m]')
grid on
% yticks([-20:10:80])


%% ===============================retreat==================================
figure('WindowState','maximized');
[AX,H1,H2] = plotyy(O.timenum(1,:),O.dndt_cliff(1,:),O.timenum(1,:),cumsum(O.dndt_cliff(1,:)),'plot');
datetick
set(AX,'FontSize',18)
set(AX(1),'FontSize',18)
set(AX(2),'FontSize',18)
set(H1,'linewidth',2)
set(H2,'linewidth',2)
title('Cliff Retreat and Accummulation Over Time in Years 2022-2052 Using TWL Model')
ylabel(AX(1),'Cliff Retreat [m/yr]')
ylabel(AX(2),'Cummulative Cliff Retreat [m]')
xlabel('Time [year]')
% set(AX(1),'YTick',[.4:0.05:0.8]);
% set(AX(2),'YTick',[0:2.5:20]);


%% ============================Drift_Vol_Shore=============================
figure('WindowState','maximized');
[AX,H1,H2] = plotyy(O.timenum(1,:),O.Drift_Vol_Shore(1,:),O.timenum(1,:),cumsum(O.Drift_Vol_Shore(1,:)),'plot');
datetick
set(AX,'FontSize',18)
set(AX(1),'FontSize',18)
set(AX(2),'FontSize',18)
set(H1,'linewidth',2)
set(H2,'linewidth',2)
title('Erosion Volume to Shore Over Time in Years 2022-2052 Using TWL Model')
ylabel(AX(1),'Erosion Volume to Shore [m3/yr]')
ylabel(AX(2),'Cummulative Erosion Volume to Shore [m3]')
xlabel('Time [year]')
% set(AX(1),'YTick',[4.5:0.25:6.5]);
% set(AX(2),'YTick',[0:25:200]);


%% ========================Drift_Vol_Shore / S.d===========================
figure('WindowState','maximized');
[AX,H1,H2] = plotyy(O.timenum(1,:),O.Drift_Vol_Shore(1,:)/S.d,O.timenum(1,:),cumsum(O.Drift_Vol_Shore(1,:)/S.d),'plot');
datetick
set(AX,'FontSize',18)
set(AX(1),'FontSize',18)
set(AX(2),'FontSize',18)
set(H1,'linewidth',2)
set(H2,'linewidth',2)
title('Shoreline Change and Accummulation Over Time in Years 2022-2052 Using TWL Model')
ylabel(AX(1),'Shoreline Change "gain" [m/yr]')
ylabel(AX(2),'Cummulative Shoreline Change "gain" [m]')
xlabel('Time [year]')
% set(AX(1),'YTick',[.45:0.025:.65]);
% set(AX(2),'YTick',[0:2.5:20]);


%% ===================== Drift_Vol_Shore / (S.d+Ej) =======================
figure('WindowState','maximized');
[AX,H1,H2] = plotyy(O.timenum(1,:),O.Drift_Vol_Shore(1,:)/(S.d+3.5),O.timenum(1,:),cumsum(O.Drift_Vol_Shore(1,:)/(S.d+3.5)),'plot');
datetick
set(AX,'FontSize',18)
set(AX(1),'FontSize',18)
set(AX(2),'FontSize',18)
set(H1,'linewidth',2)
set(H2,'linewidth',2)
title('Shoreline Change and Accummulation Over Time in Years 2022-2052 Using TWL Model')
ylabel(AX(1),'Shoreline Change "gain" [m/yr]')
ylabel(AX(2),'Cummulative Shoreline Change "gain" [m]')
xlabel('Time [year]')
% set(AX(1),'YTick',[.35:0.015:.5]);
% set(AX(2),'YTick',[0:1.5:15]);


%% ===============================retreat==================================
figure('WindowState','maximized'); 
hold on
box on
plot(O.timenum(1,:), O.dndt_cliff(1,:),'linewidth',2)
plot(O.timenum(1,:), O.dndt_cliff(8,:),'linewidth',2)
title('Cliff Retreat Over Time in Years 2022-2052 Using All Models')
xlabel('Time [Years]')
ylabel('Cliff Retreat [m/yr]')
set(gca,'FontSize',18)
set(gca, 'XTickLabel', O.timenum(1,1)-O.timenum(1,end))
% set(gca, 'XTickLabel', 2022-2055)
datetick


%% =======================retreat & accummulation==========================
figure('WindowState','maximized');
[AX,H1,H2] = plotyy(O.timenum(1,:),O.dndt_cliff(1,:),O.timenum(1,:),cumsum(O.dndt_cliff(1,:)),'plot');
datetick
% set(AX(1),'YTick',[0.76:0.005:0.82]);
% set(AX(2),'YTick',[0:2.5:30]);
set(AX,'FontSize',18)
set(AX(1),'FontSize',18)
set(AX(2),'FontSize',18)
set(H1,'linewidth',2)
set(H2,'linewidth',2)
title('Cliff Retreat and Accummulation Over Time in Years 2022-2052 Using TWL Model')
ylabel(AX(1),'Cliff Retreat [m/yr]')
ylabel(AX(2),'Cummulative Cliff Retreat [m]')
xlabel('Time [year]')


%% =============================Beach Width================================
figure('WindowState','maximized');
plot(O.timenum(1,:),O.BeachWidth(1,:),'linewidth',2);
set(gca,'FontSize',18)
title('Beach Width Over Time in Years 2022-2052 Using TWL Model')
xlabel('Time [year]')
ylabel('Beach Width [m]')
datetick


%% =======================Beach or Beach Slope=========================
figure('WindowState','maximized');
plot(O.timenum(1,:),O.tan_B(1,:),'LineWidth',2);
% legend('O.tan-B')
title('Beach Slope Over Time in Years 2022-2052 Using TWL Model')
xlabel('Time [Years]')
ylabel('Beach Slope [Meters]')
set(gca,'FontSize',18)
% set(gca, 'XTickLabel', O.timenum(1,1)-O.timenum(1,end))
% xticks([2020:2:2060])
datetick


%% =================================TWL====================================
figure('WindowState','maximized');
plot(O.timenum(1,:),O.TWL(1,:),'b','LineWidth',2);
legend('O.TWL')
title('TWL Over Time in Years 2022-2052 Using TWL Model')
xlabel('Time [Years]')
ylabel('TWL [Meters]')
set(gca,'FontSize',18)
set(gca, 'XTickLabel', O.timenum(1,1)-O.timenum(1,end))
datetick


%% ===============================eta_setup================================
figure('WindowState','maximized');
plot(O.timenum(1,:),O.eta_setup(1,:),'b','LineWidth',2);
title('Setup Over Time in Years 2022-2052 Using TWL Model')
xlabel('Time [Years]')
ylabel('Setup [Meters]')
set(gca,'FontSize',18)
set(gca, 'XTickLabel', O.timenum(1,1)-O.timenum(1,end))
datetick


%% ==========================R2_percent_Stockdon===========================
figure('WindowState','maximized');
plot(O.timenum(1,:),O.R2_percent_Stockdon(1,:),'b','LineWidth',2);
title('Runup Over Time in Years 2022-2052 Using TWL Model')
xlabel('Time [Years]')
ylabel('Runup [Meters]')
set(gca,'FontSize',18)
set(gca, 'XTickLabel', O.timenum(1,1)-O.timenum(1,end))
datetick


%% =================Beach Width / Slope / Cliff Retreat====================
figure('WindowState','maximized');
% plot(O.timenum(1,:),O.dndt_cliff(1,:));
plot(O.tan_B(1,:),O.dndt_cliff(1,:),'Linewidth',2);
% plot(O.timenum(1,:),O.BeachWidth(1,:));
set(gca,'FontSize',18)
% title('Cliff Retreat Over Time in Years 2022-2052')
title('Cliff Retreat Versus Beach Slope Over Time in Years 2022-2052')
% title('Beach Width Over Time in Years 2022-2052')
ylabel('Cliff Retreat [Meters]')
% ylabel('Beach Slope [Meters]')
% ylabel('Beach Width [Meters]')
% xlabel('Time [years]')
xlabel('Beach Slope [-]')


%% ==========================Beach Width & Slope===========================
figure('WindowState','maximized');
[AX,H1, H2] = plotyy(O.timenum(1,:),O.BeachWidth(1,:),O.timenum(1,:),O.tan_B(1,:),'plot');
datetick
set(gca,'FontSize',18)
set(AX(2),'FontSize',18)
title('Beach Width and Slope Over Time in Years 2022-2052')
ylabel(AX(1),'Beach Width [m]')
ylabel(AX(2),'Beach Slope [-]')
xlabel('Time [years]')
set(H1,'linewidth',2)
set(H2,'linewidth',2)
% set(AX(1),'YTick',[40:5:70]);
% set(AX(2),'YTick',[0.05:0.005:0.08]);


%% ========================== SL & TWL & CR ===============================
figure('WindowState','maximized');
MSL=[0 0.4 0.6 1.2 1.5 1.8 2.1]; % 5cm SLR included in TWL
TWL=[4.3050 4.6811 4.8693 5.4343 5.7172 6.0002 6.2833]; 
CR=[0.6458 0.7022 0.7304 0.8152 0.8576 0.9000 0.9425];  
% [AX,H1,H2] = plotyy(MSL,TWL,MSL,CR,'plot')
set(gca,'FontSize',18)
yyaxis left
plot(MSL,TWL,'-o','linewidth',2)
title('TWL and Cliff Retreat Versus Sea Levels Using TWL Model')
ylabel('TWL [m]')
xlabel('SL [m+MSL]')
% ylim([4 8])
% xlim([0 2.1])
% xticks([0:0.1:2.1])
yyaxis right
plot(MSL,CR,'-s','linewidth',2)
set(gca,'FontSize',18)
ylabel('Cliff Retreat [m/yr]')
% ylim([0.6 1])


%% ============================ CR VS Slope ===============================
figure('WindowState','maximized');
[AX,H1,H2] = plotyy(O.timenum(1,:),O.dndt_cliff(1,:),O.timenum(1,:),cumsum(O.dndt_cliff(1,:)),'plot');
datetick
set(AX,'FontSize',18)
set(AX(1),'FontSize',18)
set(AX(2),'FontSize',18)
set(H1,'linewidth',2)
set(H2,'linewidth',2)
title('Cliff Retreat and Accummulation Over Time in Years 2023-2051 for Transect 453')
ylabel(AX(1),'Cliff Retreat [m/yr]')
ylabel(AX(2),'Cumulative Cliff Retreate [m]')
xlabel('Time [year]')
% set(AX(1),'YTick',[0.4:0.05:.8]);
% set(AX(2),'YTick',[0:2.5:20]);


%% ============================ Beach VS Slope ============================
figure('WindowState','maximized');
[AX,H1, H2] = plotyy(O.timenum(1,:),O.BeachWidth(1,:),O.timenum(1,:),O.tan_B(1,:),'plot');
datetick
set(gca,'FontSize',18)
set(AX(2),'FontSize',18)
title('Beach Width and Slope Over Time in Years 2022-2052 for Transect 453')
ylabel(AX(1),'Beach Width [m]')
ylabel(AX(2),'Beach Slope [-]')
xlabel('Time [years]')
set(H1,'linewidth',2)
set(H2,'linewidth',2)
% set(AX(1),'YTick',[45:1.5:60]);
% set(AX(2),'YTick',[0.05:0.003:0.08]);


%% ============================ cliff position ============================
figure('WindowState','maximized');
plot(O.timenum(1,:),O.cliffposition(1,:),'linewidth',2);
datetick
set(gca,'FontSize',18)
title('Cliff Position Over Time in Years 2022-2052 for Transect 453')
xlabel('Time [year]')
ylabel('Cliff Position [m]')


%% ================================= SLR ==================================
% figure('WindowState','maximized');
% plot(SLRquad(:,1),SLRquad(:,2),'b','linewidth',2);
% datetick
% set(gca,'FontSize',18)
% title('Sea Level Rise Over Time in Years 2021-2051 for Transect 453')
% ylabel('SLR [m]')
% xlabel('Time [year]')
% % yticks([.10:0.02:.30]);


%% =============================== dndt_cliff =============================
figure('WindowState','maximized');
[AX,H1,H2] = plotyy(O.timenum(1,:),O.dndt_cliff(1,:),O.timenum(1,:),(cumsum(O.dndt_cliff(1,:))),'plot');
% [AX,H1,H2] = plotyy(O.timenum(1,:),O.dndt_cliff(1,:),SLRquad(:,1),(SLRquad(:,2)),'plot')
% [AX,H1,H2] = plotyy(O.timenum(1,:),O.dndt_cliff(1,:),O.timenum(1,:),(O.SL(1,:)),'plot')
% [AX,H1,H2] = plotyy(O.timenum(1,:),O.dndt_cliff(1,:),O.timenum(1,:),(O.TWL(1,:)),'plot')
datetick
set(AX,'FontSize',18)
set(AX(1),'FontSize',18)
set(AX(2),'FontSize',18)
set(H1,'linewidth',2)
set(H2,'linewidth',2)
title('Cliff Retreat and Total Water Level Over Time in Years 2023-2051 for Transect 453')
ylabel(AX(1),'Cliff Retreat [m/yr]')
ylabel(AX(2),'Cumulative Cliff Retreate [m]')
% % ylabel(AX(2),'SLR [m]')
% % ylabel(AX(2),'Sea Level [+SLR m]')
% ylabel(AX(2),'TWL [m]')
xlabel('Time [year]')

