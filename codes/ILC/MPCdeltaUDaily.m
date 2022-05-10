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
ut_min = 53.5;
ut_max = 65;

%Parameter Indentifiction: Forecasting Model
A = RegCoefT(1:VSign,:);
B = RegCoefT(VSign+1,:);
C = RegCoefT(VSign+2,:);

%real-time data
x = f2k(Trm);


%Real-time designed data before for energy consumption calculation
Toa = f2k (AHUdataDay(:,1));
Tma = f2k (AHUdataDay(:,4));
RAflowD = CFMCoeff * AHUdataDay(:,6);
SAflowD = CFMCoeff * AHUdataDay(:,7);
WsfD = AHUdataDay(:,11);
WrfD = AHUdataDay(:,12);

NumDayMPC = NumDayPartition - p;

%Quadratic Weight Factor
coef = [0.01 1  1.1258e+03];
% Qtr = 0.5;
% Qfr = 200;
% Qfs = 200;
% R1 = 1;
% R2 = 2;
% R3 = 400;
Qtr = 1.5;
Qfr = 222; %222
Qfs = 822; %222
R1 = coef(1);
R2 = coef(2);
R3 = coef(3);

for i = 1 : 26
    Dsign = i - 1; %Select the day for MPC
%     Sol(:,Dsign*NumDayMPC+1:Dsign*NumDayMPC+2)= [Tra(Dsign*NumDayMPC+1:Dsign*NumDayMPC+2) RAflow(Dsign*NumDayMPC+1:Dsign*NumDayMPC+2)...
%         betaD(Dsign*NumDayMPC+1:Dsign*NumDayMPC+2) Tsa(Dsign*NumDayMPC+1:Dsign*NumDayMPC+2) ...
%         OAflow(Dsign*NumDayMPC+1:Dsign*NumDayMPC+2)]';
    for j = 1 : NumDayPartition - p - 2
        Hsign = j ; %Select the hour after occupied mode turns on(7 am). Ex: 2 * 4 means 9 am for 4 samples 1 hour

        %Define the disigned data in MPC steps
        % Ddata = [Toa RAflowD SAflowD WsfD WrfD];
        % D = Ddata(Dsign * NumDayMPC + Hsign : Dsign * NumDayMPC + Hsign + 1, :);
        %Setpoint of room temperature and SA flow
        WPrev = Toa(Dsign * NumDayMPC + Hsign);
        WCurr = Toa(Dsign * NumDayMPC + Hsign + 1);
        WPred = Toa(Dsign * NumDayMPC + Hsign + 2);
        RSratio = RAflowD(Dsign * NumDayMPC + Hsign + 2) / SAflowD(Dsign * NumDayMPC + Hsign + 1);
        SP = [TrmReg(Dsign * NumDayMPC + Hsign + 2)
            SAflowIn(Dsign * NumDayMPC + Hsign + 1)
            RSratio * SAflowIn(Dsign * NumDayMPC + Hsign + 1)];

        XtPrev = xt0(Dsign * NumDayMPC + Hsign,:);
        XfPrev = xf0(Dsign * NumDayMPC + Hsign,:);

        UPrev = u0(Dsign * NumDayMPC + Hsign,:)';

        XtCurr = aT * XtPrev + bT * UPrev + cT * WPrev;
        XfCurr = aF * XfPrev + bF * UPrev + cF * WPrev;

        %Jocobian
        deltaJ = [ [ [Qtr 0;0 Qfr] zeros(2,3) eye(2)];
            zeros(3,2) [-2*R1-Qfs*2*XfCurr^2 Cv * XfCurr -Qfs*2*XfCurr; Cv * XfCurr -2*R2 Cv; -Qfs*2*XfCurr Cv -Qfs*2-2*R3] [bT' bF'];
            [eye(2) [-bT;-bF] zeros(2)] ];


        B = [Qtr*SP(1) Qfr*SP(3) Cv*XfCurr*XtCurr+2*Qfs*XfCurr*(UPrev(1)*XfCurr+UPrev(3)-SP(2)) -Cv*(UPrev(3)+UPrev(1)*XfCurr)...
            Cv*(WCurr-UPrev(2))+2*Qfs*(UPrev(1)*XfCurr+UPrev(3)-SP(2)) aT*XtCurr+bT*UPrev+cT*WCurr aF*XfCurr+bF*UPrev+cF*WCurr]';

        DV = deltaJ \ B;

        DV(3:5) = DV(3:5) + UPrev;

        % Saturation
        if(DV(3)>u1max)
            DV(3)=u1max;
        end
        if(DV(4)>u2max)
            DV(4)=u2max;
        end
        if(DV(5)>u3max)
            DV(5)=u3max;
        end
        if(DV(3)<u1min)
            DV(3)=u1min;
        end
        if(DV(4)<u2min)
            DV(4)=u2min;
        end
        if(DV(5)<u3min)
            DV(5)=u3min;
        end

        Sol( : , (NumDayMPC - 2) * Dsign + Hsign) = DV;
    end
end

save MPCresultDataDaily.mat

AHU9MPCNumericalMayProcessPlotDaily
