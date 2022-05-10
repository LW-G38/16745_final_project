clear all
close all

load VAVMPCDailyForecast.mat
% back to full length data
x0 = Trm;
u_ref = Tdisch.*flow;
w_forecast = Toa_forecast;
w = Toa;
tt = linspace(0, DayNum, length(xt));
N = length(tt);

%% torque limits
% airflow limits
uf_2f_min = [200 190 160 110 120 130 110 110 120 130 110 75 60 110 110 510];
uf_2f_max = [800 750 630 450 470 530 450 450 470 530 450 300 230 440 430 2030];

uf_3f_min = [195 170 150 110 120 120 110 110 120 120 110 75 60 105 105 225];
uf_3f_max = [780 690 600 450 470 490 450 450 470 490 450 300 230 420 420 1020];

% discharge air temperature limits
ut_min = 55;
ut_max = 70;

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
Tsp = f2k(VData(:,(1:32) + 2*32));
%%
p = 3;
NumDayMPC = NumDayPartition - p;
%Quadratic Weight Factor
Q = 1000000*eye(32);
R = 0.01*eye(32);

%accessory
I = eye(32);
O = zeros(32,32);

% Hqp = blkdiag(Q,R);
Hqp = blkdiag(Q,Q,Q,R,R,R);
for i = 1 : 26
    Dsign = i - 1; %Select the day for MPC
    for j = 1 : NumDayMPC
        Hsign = j ; %Select the hour after occupied mode turns on(7 am). Ex: 2 * 4 means 9 am for 4 samples 1 hour

        %Define the disigned data in MPC steps
%         WPrev = Toa_forecast(Dsign * NumDayMPC + Hsign);
%         XtPrev = xt0(Dsign * NumDayMPC + Hsign,:)';
%         UPrev = u0(Dsign * NumDayMPC + Hsign,:)';
%         Xref = Tsp(Dsign * NumDayMPC + Hsign + 1,:)';
        Wf = Toa_forecast(Dsign * NumDayMPC + Hsign: Dsign * NumDayMPC + Hsign + 2);
        Xprev = xt0(Dsign * NumDayMPC + Hsign,:)';
        Xref = Tsp(Dsign * NumDayMPC + Hsign + 1 : Dsign * NumDayMPC + Hsign + 3,:)';        
%         Fqp = [-Q'*Xref; zeros(32,1)];
%         Aeq = [eye(32) -B];
%         beq = [A*XtPrev + C*WPrev];
%         lb = [-Inf*ones(32,1); u_min'];
%         ub = [Inf*ones(32,1); u_max'];
        Fqp = [-Q'*Xref(:,1); -Q'*Xref(:,2); -Q'*Xref(:,3); zeros(32*3,1)];
        Aeq = [-I O O B O O; A -I O O B O; O A -I O O B];
        beq = -[A*Xprev + C*Wf(1); C*Wf(2); C*Wf(3)];
        lb = [-Inf*ones(32*3,1); u_min'; u_min'; u_min'];
        ub = [Inf*ones(32*3,1); u_max'; u_max'; u_max'];
        Sol(:, NumDayMPC * Dsign + Hsign) = quadprog(Hqp,Fqp,[],[],Aeq,beq,lb,ub);
    end
end

%%
Xmpc = Sol(1:96, :);
%save MPCresultDataDaily.mat