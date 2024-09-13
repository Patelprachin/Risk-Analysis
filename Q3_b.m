clear all
close all
clc
format short

rng(123)

% Parameters
FV = 100;
CR = 0.03;
annual_coupon = FV * CR;
maturity = 5;
P0 = 99;
num_simulations = 100000;
num_days = [1 10 20 30 40 50 60 70 80 90];
daily_std_dev = 0.006;
daily_mean = 0;

% Calculate coupon payments
coupon_payments = annual_coupon * ones(1, maturity);

% Calculating YTM
syms x
yield = vpasolve(annual_coupon*sum(1./(1+x).^[1:5])+FV/(1+x)^5==P0,x);
y = double(round(yield(1,1),5));

% Calculating zscore
z_score = norminv(0.01,0,1);

% Calculating VaR - Exact Formula
for i = 1:length(num_days)
    VaR_exact(i) = z_score*daily_std_dev*sqrt(num_days(1,i))*P0;
end

% Calculating ES - Exact Formula
for i = 1:length(num_days)
    ES_exact(i) = -normpdf(z_score)*daily_std_dev*sqrt(num_days(1,i))*P0/0.01;
end

% Combine VaR and ES results into a table
VaR_methods = {'Exact Formula'};
ES_methods = {'Exact Formula'};
VaR_results = VaR_exact';
ES_results = ES_exact';
table_results = array2table([num_days' VaR_results, ES_results], 'VariableNames', {'Days','VaR', 'ES'});

% Add number of days as row names
% table_results.Properties.RowNames = arrayfun(@num2str, num_days, 'UniformOutput', false);

% Display the table
disp('Table comparing VaR and ES results for different number of days:');
disp(table_results);