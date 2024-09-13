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
disp('Yield to Maturity =') 
fprintf('%.4f\n', yield(1,1));

% Calculating YTM for 10% decline in bond over 30 days
syms x
yield_30 = double(vpasolve(annual_coupon*sum(1./(1+x).^[1:5])+FV/(1+x)^5==P0*0.9,x));

% For bond price to fall below 10% over 30 days yield should be greater than equal to yield_30
prob1 = normcdf(yield_30(1,1),daily_mean,daily_std_dev*sqrt(30));

% Probability of 10% decline in 30 days
prob = 1 - prob1;
disp('Probability of 10% decline in 30 days')
fprintf('%.4f\n', prob);

% Calculating zscore
z_score = norminv(0.01,0,1);

% Calculating VaR - Exact Formula
for i = 1:length(num_days)
    VaR_exact(i) = z_score*daily_std_dev*sqrt(num_days(1,i))*P0;
end

% Calculating VaR - delta approximation exact formula
y = double(round(yield(1,1),5));
dur = (annual_coupon*sum([1:5].*1./(1+y).^[1:5])+ (5*FV)/(1+y)^5)/P0;
for i = 1:length(num_days)
    deltay(i) = randn(1,1)*0.006*sqrt(num_days(i)) + 0;
    P1(i) = P0 + P0*dur*deltay(i)/(1+y);
    VaR_delta(i) = z_score*daily_std_dev*sqrt(num_days(1,i))*P1(i);
end

% Calculating VaR - delta-gamma approximation exact formula
conv = (annual_coupon*sum([1:5].*[2:6].*1./(1+y).^[1:5])+ (5*6*FV)/(1+y)^5)/P0*(1+y)^2;
for i = 1:length(num_days)
    P2(i) = P0 + P0*dur*deltay(i)/(1+y) + P0*0.5*conv*deltay(i)^2;
    VaR_gamma(i) = z_score*daily_std_dev*sqrt(num_days(1,i))*P2(i);
end

% Simulate variations in yield to maturity
ytm_simulations = zeros(length(num_days), num_simulations);
for h = 1:length(num_days)
    for i = 1:num_simulations
        ytm_simulations(h,i) = randn() * daily_std_dev * sqrt(num_days(h)) + daily_mean;
    end
end

ytm_sim = y + ytm_simulations;
ytm = max(ytm_sim, 0);

% Calculate bond prices using simulated yield to maturities
bond_prices = zeros(length(num_days), num_simulations);
for h = 1:length(num_days)
    for i = 1:num_simulations
    bond_prices(h,i) = sum(coupon_payments ./(1+ytm(h,i)).^(1:1:5)) + FV / (1+ytm(h,i)).^maturity;
    end
end

% VaR - Monte Carlo Full Revaluation
for j = (1:1:10)
    % Full Revaluation
    for f = 1:num_simulations
        fullreval_VaR(j,f) = bond_prices(j,f)*norminv(0.01)*sqrt(num_days(1,j))*daily_std_dev;
    end
end

results = sum(fullreval_VaR,2)/num_simulations;

% VaR - Monte Carlo Delta Approximation
for h = 1:length(num_days)
    for i = 1:num_simulations
        dy(h,i) = ytm(h,i) - y;
        bond_delta(h,i) = P0 + P0*dur*dy(h,i)/(1+y);
        VaR_delta_mc(h,i) = bond_delta(h,i)*norminv(0.01)*sqrt(num_days(1,h))*daily_std_dev;
    end
end

res_delta = sum(VaR_delta_mc,2)/num_simulations;

% VaR - Monte Carlo Delta-Gamma Approximation
for h = 1:length(num_days)
    for i = 1:num_simulations
        dy(h,i) = ytm(h,i) - y;
        bond_gamma(h,i) = P0 + P0*dur*dy(h,i)/(1+y) + P0*0.5*conv*dy(h,i)^2;
        VaR_gamma_mc(h,i) = bond_gamma(h,i)*norminv(0.01)*sqrt(num_days(1,h))*daily_std_dev;
    end
end

res_gam = sum(VaR_gamma_mc,2)/num_simulations;

% Create table for VaR comparison
VarTable = table(num_days', VaR_exact', VaR_delta', VaR_gamma', 'VariableNames', {'Days', 'Exact Formula', 'Delta Approximation', 'Delta Gamma Approximation'});
disp(VarTable);

% Create table for Monte Carlo comparison
MonteCarloTable = table(num_days', results, res_delta, res_gam, 'VariableNames', {'Days', 'Full Revaluation MC', 'Delta Approximation MC', 'Delta Gamma Approximation MC'});
disp(MonteCarloTable);