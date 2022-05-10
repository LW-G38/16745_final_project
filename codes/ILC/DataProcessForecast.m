clear all
close all
% load VAVdataMay2014.mat % vav data lacks 5/13, 5/17, 5/18, 5/19
filename = '118030_May_2014.csv';

table_forecast = readtable(filename);

%%
oat_forecast = table_forecast.Temperature;

year_data = table_forecast.Year;
month_data = table_forecast.Month;
date_data = table_forecast.Day;
% hour_data = table_forecast.Hour;
% minute_data = table_forecast.Minute;

date_vec = [year_data month_data date_data];

t = datetime(date_vec);

oat_forecast(t == ("13-May-2014") | t== ("17-May-2014") | t== ("18-May-2014") | t== ("19-May-2014")) = [];
% t(t == ("13-May-2014") | t== ("17-May-2014") | t== ("18-May-2014") | t== ("19-May-2014")) = [];
oat_forecast = repelem(oat_forecast,2);

clear t table_forecast date_vec date_data year_data month_data filename
save forecast.mat