clear all
close all
swEPSfigure
load MPCresultDataDaily.mat

TrmCalc = Sol(1:32,:);
uMPC = Sol(33:64,:);

%Parameter Indentifiction: Forecasting Model
A = RegCoefT(1:VSign,:)';
B = RegCoefT(VSign+1,:)';
B = diag(B);
C = RegCoefT(VSign+2,:)';

Atrue = RegCoefW(1:VSign,:)';
Btrue = diag(RegCoefW(VSign+1,:)');
Ctrue = RegCoefW(VSign+2,:)';

TrmFore = zeros(size(TrmCalc));
xt0 = xt0';

for i = 1 : DayNum
    Dsign = i - 1;
    TrmFore(:, Dsign * NumDayMPC+1) = xt0(:,Dsign * NumDayMPC+1);
    for j = 2 : NumDayMPC
        TrmPrev = xt0(:,Dsign * NumDayMPC+j-1);
        TrmFore(:, Dsign * NumDayMPC + j) = Atrue*TrmPrev+...
            + Btrue*uMPC(:,Dsign * NumDayMPC+j-1) + Ctrue*w(Dsign * NumDayMPC+j-1);
    end
end



%%
%Discard every first two samples of measurements
NumDayMPC  = NumDayMPC - 2;
for i = 1 : DayNum
    OAflow(NumDayMPC * (i - 1) + 1 : NumDayMPC * (i - 1) + 2) = [];
    betaD(NumDayMPC * (i - 1) + 1 : NumDayMPC * (i - 1) + 2) = [];
    RAflow(NumDayMPC * (i - 1) + 1 : NumDayMPC * (i - 1) + 2) = [];
    SAflowD(NumDayMPC * (i - 1) + 1 : NumDayMPC * (i - 1) + 2) = [];
    Tma(NumDayMPC * (i - 1) + 1 : NumDayMPC * (i - 1) + 2) = [];
    Tsa(NumDayMPC * (i - 1) + 1 : NumDayMPC * (i - 1) + 2) = [];
    Tra(NumDayMPC * (i - 1) + 1 : NumDayMPC * (i - 1) + 2) = [];
    w(NumDayMPC * (i - 1) + 1 : NumDayMPC * (i - 1) + 2) = [];
    TrmReg(NumDayMPC * (i - 1) + 1 : NumDayMPC * (i - 1) + 2) = [];
    SAflowIn(NumDayMPC * (i - 1) + 1 : NumDayMPC * (i - 1) + 2) = [];
    RAflowD(NumDayMPC * (i - 1) + 1 : NumDayMPC * (i - 1) + 2) = [];
    Toa(NumDayMPC * (i - 1) + 1 : NumDayMPC * (i - 1) + 2) = [];
end


En = Cv * OAflow .* w + Cv * betaD .* RAflow .* Tra - Cv * (SAflowD) .* Tsa;
En = Cv * SAflowD .* (Tma - Tsa);
Enopt = Cv * u3 .* w' + Cv * u1 .* x1 .* x2 - Cv * (u3 + u1 .* x2) .* u2;
Enopt = Enopt';

for i = 1:length(En)
if(En(i)<0 && Enopt(i)<0)
   En(i)=0.01; 
   Enopt(i)=0.01;
end
if(Enopt(i)<0 && En(i)>0)
    Enopt(i)=En(i); 
end
if(Enopt(i)>0 && En(i)<0)
    En(i)=Enopt(i); 
end
end
EnSave = En - Enopt;

EnPct = EnSave ./ En * 100;
for i = 1:length(En)
if(abs(EnPct(i))>=100)
    EnPct(i)=100;
end
end

for i = 1 : DayNum
    EnPctDaily(i) = mean( EnPct((i-1) * NumDayMPC + 1 : i * NumDayMPC) );
    EnSaveDaily(i) = mean( EnSave((i-1) * NumDayMPC + 1 : i * NumDayMPC) );
    EnDaily(i) = mean( En((i-1) * NumDayMPC + 1 : i * NumDayMPC) );
    EnoptDaily(i) = mean( Enopt((i-1) * NumDayMPC + 1 : i * NumDayMPC) );
end
EnPctDaily = EnPctDaily';
EnSaveDaily = EnSaveDaily';
EnDaily = EnDaily';
EnoptDaily = EnoptDaily';
% figure
% bar(EnPctDaily)
% xlim([0 27])

IndexPeak = find(En>=40);
EnoptPeak = Enopt(IndexPeak);
EnPeak = En(IndexPeak);
EnPctPeak = EnPct(IndexPeak);

x1 = k2f(x1);
x2 = x2 / CFMCoeff;
u2 = k2f(u2);
u3 = u3 / CFMCoeff;
u1Peak = u1(IndexPeak);
x2Peak = x2(IndexPeak);
u3Peak = u3(IndexPeak);

u1minhist = u1min * ones(DayNum * NumDayMPC,1);
u1maxhist = u1max * ones(DayNum * NumDayMPC,1);
u2minhist = k2f(u2min) * ones(DayNum * NumDayMPC,1);
u2maxhist = k2f(u2max) * ones(DayNum * NumDayMPC,1);
u3minhist = u3min / CFMCoeff * ones(DayNum * NumDayMPC,1);
u3maxhist = 8000 * ones(DayNum * NumDayMPC,1);


t = linspace(0,26,length(En));
% figure
% plot(t,k2f(Toa))
% xlabel('time')
% ylabel('Fresh Air Temp /F')

figure
subplot(3,1,1)
plot(t,u1maxhist,'g-')
hold on
p1=plot(t,u1);
xlim([0 26])
ylim([0 1])
hold on
p2=plot(t,betaD,'r--');
hold on
plot(t,u1minhist,'g-')
ylabel('RA ratio')
legend([p1,p2],'MPC','Original','Location','SouthEast')
subplot(3,1,2)
plot(t,u2maxhist,'g-')
hold on
plot(t,u2)
xlim([0 26])
ylim([50 75])
hold on
plot(t,k2f(Tsa),'--r')
hold on
plot(t,u2minhist,'g-')
ylabel('$T_{SA}$ (F)')
subplot(3,1,3)
plot(t,u3maxhist,'g-')
hold on
plot(t,u3)
xlim([0 26])
ylim([1800 10000])
hold on
plot(t,OAflow / CFMCoeff,'--r')
hold on
plot(t,u3minhist,'g-')
xlabel('Days (May)')
ylabel('$\dot{m_{OA}}$ (cfm)')
print -depsc FigControlMay.eps

figure
subplot(3,1,1)
plot(t,k2f(TrmReg),'r-.')
hold on
plot(t,x1)
xlim([0 26])
ylim([70 78])
ylabel('$T_{RA}$ (F)')
legend('MPC','Original','Location','NorthEast')
subplot(3,1,2)
plot(t,SAflowIn / CFMCoeff,'r-.')
hold on
plot(t,u1.*x2+u3)
xlim([0 26])
ylim([3000 15000])
ylabel('$\dot{m_{SA}}$ (CFM)')
subplot(3,1,3)
plot(t,RAflowD / CFMCoeff,'r-.')
hold on
plot(t,x2)
xlim([0 26])
ylim([0 12000])
xlabel('Days (May)')
ylabel('$\dot{m_{RA}}$ (CFM)')
print -depsc FigstatedailyMay.eps

figure
scatter(u1.*x2, u3, 20, EnPct ,'s', 'filled')
xlabel('Recirculation Air flow rate / CFM')
ylabel('Outside Air flow rate / CFM')
caxis([-20, 100])
colorbar

figure
scatter(u1.*x2, u3, 20, Enopt/max(En) ,'s', 'filled')
xlabel('Recirculation Air flow rate / CFM')
ylabel('Outside Air flow rate / CFM')
ylim([1500 8000])
caxis([0, 1])
colorbar
print -depsc FigEnoptScatter.eps

figure
scatter(betaD.*RAflowD / CFMCoeff, OAflow / CFMCoeff, 20, En/max(En) ,'s', 'filled')
xlabel('Recirculation Air flow rate / CFM')
ylabel('Outside Air flow rate / CFM')
ylim([1500 8000])
caxis([0, 1])
colorbar
print -depsc FigEnScatter.eps

figure
subplot(212)
scatter(u1.*x2, u3, 20, Enopt/max(En) ,'s', 'filled')
title('MPC')
ylabel({'Outside Air';' Flow Rate (cfm)'})
xlabel('Recirculation Air Flow Rate (cfm)')
ylim([1500 8000])
xlim([1500 10000])
caxis([0, 1])
subplot(211)
scatter(betaD.*RAflowD / CFMCoeff, OAflow / CFMCoeff, 20, En/max(En) ,'s', 'filled')
title('Original')

ylabel({'Outside Air';' Flow Rate (cfm)'})
ylim([1500 8000])
xlim([1500 10000])
caxis([0, 1])
print -depsc FigEnoptEnScatter.eps

tday = linspace(7.75, 25, NumDayPartition - p - 2);
for i = 1 : DayNum
figure
plot(tday,En((i-1) * ( NumDayPartition - p - 2) + 1 : i * ( NumDayPartition - p - 2)),'r*-')
hold on
plot(tday,Enopt((i-1) * ( NumDayPartition - p - 2) + 1 : i * ( NumDayPartition - p - 2)))
xlabel('time')
ylabel('Energy/kJ')
%print -depsc energy_trend.eps
end

figure

plot(t,En,'r-.',t,Enopt)
xlim([0 DayNum])
xlabel('Days (May)')
ylabel('Energy/kJ')
legend('Original','MPC')
print -depsc energy_trend.eps



%Day to Day Analysis
%Daily Parameter: OAflow, RAflow, SAflow, Toa, Tra, Tsa, Tra from VAVs by Linear Regression, betaD
%Optimized Tra (x1), Optimized RAflow(x2),Optimized Tsa (u2), Optimized OAflow(u3),
%Optimized SAflow(u1.*x2+u3), Optimized beta
%Difference between Tra, RAflow, Tsa, OAflow, SAflow, beta
Para = [OAflow RAflow SAflowD Toa Tra Tsa TrmReg betaD x1' x2' u2' u3' (u1.*x2+u3)' u1' ...
    x1'-Tra x2'-RAflow u2'-Tsa u3'-OAflow (u1.*x2+u3)'-SAflowD u1'-betaD];

for i = 1 : DayNum
    ParaDaily(i,:) = mean( Para((i - 1) * NumDayMPC + 1 : i * NumDayMPC, : ) );
end

%Correlation from Parameters to Energy Savings
CorrEnPct = corr(EnPctDaily, ParaDaily);
CorrEnSave = corr(EnSaveDaily, ParaDaily);
CorrEn = corr(EnDaily, ParaDaily);
CorrEnopt = corr(EnoptDaily, ParaDaily);

CorrTable = [CorrEn; CorrEnopt; CorrEnSave; CorrEnPct];

%To daily Toa
ToaDaily = ParaDaily(:,4);
figure
[EnToa,EnToaLine1,EnToaLine2] = plotyy(1:DayNum,EnPctDaily,1:DayNum,k2f(ToaDaily),'bar','plot');

set(EnToaLine1,'FaceColor','y','EdgeColor','k');
set(EnToaLine2,'Marker','s','MarkerEdgeColor','k','MarkerFaceColor','g');

xlabel('Days (May)')
ylabel(EnToa(1),'Energy Saving $\%$') 
ylabel(EnToa(2),'$T_{OA}$ (F)') 
xlim(EnToa(1),[0 27])
xlim(EnToa(2),[0 27])

legend('Energy Saving $\%$','$T_{OA}$ (F)')

print -depsc EnPctToa.eps

RAflowDaily = ParaDaily(:,2);
figure
[EnRAflow,EnRAflowLine1,EnRAflowLine2] = plotyy(1:DayNum,EnPctDaily,1:DayNum,RAflowDaily / CFMCoeff,'bar','plot');

set(EnRAflowLine1,'FaceColor','y','EdgeColor','k');
set(EnRAflowLine2,'Marker','s','MarkerEdgeColor','k','MarkerFaceColor','g');

xlabel('Days (May)')
ylabel(EnRAflow(1),'Energy Saving $\%$') 
ylabel(EnRAflow(2),'$\dot{m_{RA}}$ / CFM') 
xlim(EnRAflow(1),[0 27])
xlim(EnRAflow(2),[0 27])

legend('Energy Saving $\%$','$\dot{m_{RA}}$ / CFM')

print -depsc EnergyRAflow.eps

%xlswrite('CorrelationMatrix.xlsx', CorrTable,'A1:T4')

%save MayFinal.mat

hourspan = 8: 1 : 25;
timespan = [hourspan(1)];
for h = 2 : length(hourspan)
timespan = [timespan hourspan(h) * ones(1,4)];
end
Thour = repmat(timespan,1,DayNum) ;

figure
scatter(u1.*x2, u3, 20, Thour ,'s', 'filled')
xlabel('Recirculation Air Flow rate / CFM')
ylabel('Outside Air flow rate / CFM')
ylim([1500 8000])
caxis([8, 25])
colorbar

figure
scatter(betaD.*RAflowD / CFMCoeff, OAflow / CFMCoeff, 20, Thour ,'s', 'filled')
xlabel('Recirculation Air flow rate / CFM')
ylabel('Outside Air flow rate / CFM')
ylim([1500 8000])
caxis([8, 25])
colorbar

figure
scatter3(u1.*x2, u3, Thour,20,En'/max(En),'s', 'filled')
figure
scatter3(betaD.*RAflowD / CFMCoeff, OAflow / CFMCoeff, Thour' ,20,En/max(En),'s', 'filled')

