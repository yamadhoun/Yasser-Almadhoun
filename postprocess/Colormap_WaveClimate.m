global S
tic
Hso=0:0.4:4;
tper=linspace(0,15,size(Hso,2));
m=length(Hso);
n=length(tper);
for i=1:m
    S.Hso = Hso(i);
    for j=1:n
        S.tper = tper(j);
        run ShorelineS_test1.m
        copy = cumsum(O.dndt_cliff(11,:));
        Ozz(i,j) = copy(end)
    end
end


toc

%% plot 1
figure; plot((Ozz)); colorbar
title('Cliff Retreat Verus Wave Climate Over Time in Years 2022-2052 Using TWL Model')
xlabel('time step [yr]')
ylabel('cliff retreat [m/yr]')
set(gca,'FontSize',18)

%% plot 2
% figure; contour(Hso,tper,squeeze(Ozz)); colorbar
figure; box on; surf(Hso,tper,squeeze(Ozz)); view(0,90)
hcb=colorbar; % hcb.Title.String = "Cummulative Cliff Retreat after 31 years [m]";
ylabel(hcb, 'Cummulative Cliff Retreat after 31 years [m]')
% set(hcb,'XTick',[5:2.5:25])
title('Cliff Retreat Verus Wave Climate (Hs,Tp) Over Time in Years 2022-2052 Using TWL Model')
xlabel('Significant Wave Height Hs [m]')
ylabel('Peak Period Tp [s]')
set(gca,'FontSize',18)
% xticks([0:0.25:4])
% yticks([0:1.5:15])
shading interp




toc