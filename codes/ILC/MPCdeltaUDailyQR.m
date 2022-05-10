clear all
close all

load VAVMPCDailyForecast.mat

%% torque limits
% airflow limits
uf_2f_min = [200 190 160 110 120 130 110 110 120 130 110 75 60 110 110 510];
uf_2f_max = [800 750 630 450 470 530 450 450 470 530 450 300 230 440 430 2030];

uf_3f_min = [195 170 150 110 120 120 110 110 120 120 110 75 60 105 105 225];
uf_3f_max = [780 690 600 450 470 490 450 450 470 490 450 300 230 420 420 1020];

% discharge air temperature limits
ut_min = 50;
ut_max = 80;

% control limits

u_min = f2k(ut_min)*CFMCoeff*[uf_2f_min uf_3f_min];
u_max = f2k(ut_max)*CFMCoeff*[uf_2f_max uf_3f_max];

%Parameter Indentifiction: Forecasting Model
A = RegCoefT(1:VSign,:)';
B = RegCoefT(VSign+1,:)';
B = diag(B);
C = RegCoefT(VSign+2,:)';

%Real-time designed data before for energy consumption calculation
Toa_forecast = w_forecast;
Tsp = f2k(VData(:,3*(1:32)));
%%
NumDayMPC = NumDayPartition - p;

%Quadratic Weight Factor
Q = eye(32);
R = 0.01*eye(32);

for i = 1 : 26
    Dsign = i - 1; %Select the day for MPC
    for j = 1 : NumDayPartition - p - 2
        Hsign = j ; %Select the hour after occupied mode turns on(7 am). Ex: 2 * 4 means 9 am for 4 samples 1 hour

        WPrev = Toa_forecast(Dsign * NumDayMPC + Hsign);
        WCurr = Toa_forecast(Dsign * NumDayMPC + Hsign + 1);

        XtPrev = xt0(Dsign * NumDayMPC + Hsign,:)';
        UPrev = u0(Dsign * NumDayMPC + Hsign,:)';

        XtCurr = A*XtPrev + B*UPrev + C*WPrev;
        Xref = Tsp(Dsign * NumDayMPC + Hsign + 2,:)';

        %Jocobian
        deltaJ = [Q*Xref; -R*UPrev; A*XtCurr+B*UPrev+C*WCurr];

        MV = [Q zeros(32,32) eye(32); zeros(32,32) R B'; eye(32) -B zeros(32,32)];
        
        DV = deltaJ \ MV;

        Sol(:, (NumDayMPC - 2) * Dsign + Hsign) = DV;

        Ucalc = DV(33:64) + UPrev';

        Ucalc = clamp(Ucalc, u_min, u_max);

        Sol(33:64, (NumDayMPC - 2) * Dsign + Hsign) = Ucalc;
    end
end

save MPCresultDataDaily.mat