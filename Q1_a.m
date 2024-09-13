clear all
close all
clc
format short
%**************************************************************************
%% *GET DATA * 
%**************************************************************************
Stock_names = ["INTC","JPM","AA","PG","MSFT"];
start_date = '01-Jan-2014';
end_date = '31-Dec-2023';
Dates = get_MarketDataViaYahoo('INTC',start_date, end_date,'1d');
df = [Dates(:,1) ];
for i = Stock_names
    df_values = get_MarketDataViaYahoo(convertStringsToChars(i),start_date, end_date,'1d');
    df_adj = df_values(:,6); %Index for col = 6 because it's the 6TH scraped column from YFinance
    df_adj.Properties.VariableNames = i;
    df = [df df_adj]; 
end
%***************************************************************************
%% *COMPUTE LOG-RETURNS* 
%***************************************************************************
prices = table2array(df(:,2:end));
ret = zeros(height(prices),size(prices,2));
for i = 2:size(prices,1)
    for j = 1:size(prices,2)
        ret_val(i,j) = log((prices(i,j))./ (prices(i-1,j)));
        ret = ret_val;
        
    end

    
end
ret(1,:)= []; %drop the first row if you want jejeje
%***************************************************************************
%% *BUILD AN EQUALLY WEIGHTED PORTFOLIO* 
%***************************************************************************
weight = 1 / size(ret,2);
weights = repmat(weight, 1, size(ret, 2));
EW_Port = sum(ret.* weights,2);                %check on excel as well as ret

for j = 1:size(ret,2)
    figure;
    plot(ret(:,j));
    ylabel(Stock_names(j))
    xlabel('Time')
end
%***************************************************************************
%% *DESCRIPTIVE STATISTICS* 
%***************************************************************************
data = EW_Port; %Modify right hand side of the equation for Summary Stats of different df
Descriptive = struct(...
    'Number_of_Observations', size(data, 1), ...
    'Max', max(data(:)), ...
    'Min', min(data(:)), ...
    'Median', median(data(:)), ...
    'Mean', mean(data(:)), ...
    'Variance', var(data(:)), ...
    'Standard_Deviation', std(data(:)),...
    'Skewness', skewness(data(:)), ... 
    'Kurtosis', kurtosis(data(:))...
);
Descriptive_table = struct2table(Descriptive)
outside = Descriptive_table.Mean+2*Descriptive_table.Standard_Deviation;
less_than=100*((sum(data<-outside)+sum(data>outside))/Descriptive_table.Number_of_Observations)
%*************************************************************************
%% *THEORETICAL VS OBSERVED* 
%*************************************************************************
% Calculate histogram parameters
numBins = round(sqrt(Descriptive.Number_of_Observations));
histEdges = linspace(min(data), max(data), numBins);

% Gaussian PDF
x = linspace(min(data), max(data), 1000);
gaussianPDF = normpdf(x, Descriptive.Mean, Descriptive.Standard_Deviation);

% Plot
figure('Color',[1 1 1]);

% Plot main histogram

histogram(data, histEdges, 'Normalization', 'pdf')
hold on
plot(x, gaussianPDF, 'r', 'LineWidth', 2)
xlabel('Returns','interpreter','latex')
title('Equally Weighted Portfolio Returns : Empirical vs Gaussian Density','interpreter','latex')
xlim([Descriptive.Mean-6*Descriptive.Standard_Deviation Descriptive.Mean+6*Descriptive.Standard_Deviation])

figure('Color',[1 1 1]);
subplot(1,2,1)
histogram(data, histEdges, 'Normalization','pdf')
hold on
plot(x, gaussianPDF, 'r', 'LineWidth', 2)
leftTailRange = [prctile(data,1), prctile(data,45)];
xlim(leftTailRange)
xlabel('Returns','interpreter','latex')
title('Equally Weighted Portfolio Returns Left Tail','interpreter','latex')

% Plot right tail
subplot(1,2,2)
histogram(data, histEdges, 'Normalization','pdf')
hold on
plot(x, gaussianPDF, 'r', 'LineWidth', 2)
rightTailRange = [prctile(data,55), prctile(data,99)];
xlim(rightTailRange)
xlabel('Returns','interpreter','latex')
title('Equally Weighted Portfolio Returns Right Tail','interpreter','latex')

%*************************************************************************
%% *QQ-PLOT*
%*************************************************************************
data_standardized = ((data-Descriptive.Mean)./Descriptive.Standard_Deviation);
quant = lognrnd(Descriptive.Mean, Descriptive.Standard_Deviation, 2515, 1);
figure
qqplot(data)
title('QQ-Plot', 'interpreter','latex')
xlabel('Standard Normal Quantiles','interpreter','latex')
ylabel('Portfolio Quantiles','interpreter','latex')
%*************************************************************************
%% *JARQUE - BERA TEST* 
%*************************************************************************
JB1 = (Descriptive.Number_of_Observations/6.*Descriptive.Skewness^2)+...
    (Descriptive.Number_of_Observations/24.)*((Descriptive.Kurtosis-3)^2);
data_standardized = (data-Descriptive.Mean)./Descriptive.Standard_Deviation;
JB2 = jbtest(data_standardized,0.05); % IF JB2 equal to 1, reject null hypothesis
p_JB = chi2cdf(JB1,2)
%Above repoirted the first method f evaluating the JB test since it's not
%good for small sample sizes, we provide below the matlab corrected for N
%of observation JB test
%To find chisquare Critical value with 2 degrees of freedom and 0.95 level
%of significance use theinverse cumulative function such that:
CV = chi2inv( 0.95, 2);

if JB1 > CV
    disp('JB test with original formulation states that data is not Normal')
else
    disp('JB test with original formulation states that data is Normal')
end

if JB2 == 1
    disp('JB test with formulation with correction for Number of Observations states that data is not Normal')
else
    disp('JB test with formulation with correction for Number of Observations states that data is Normal')
end
Jarque_Bera_table = struct(...
    'Jarque_Bera_test', JB1,...
    'pValue',chi2pdf(JB1,2));
Jarque_Bera_table = struct2table(Jarque_Bera_table)
%%
window = 120;
n_iterations = floor(size(data,1) / window); % Numero di iterazioni senza sovrapposizione
h = zeros(n_iterations, 1);
p = zeros(n_iterations, 1);
jbstat = zeros(n_iterations, 1);
critval = zeros(n_iterations, 1);

for i = 1:n_iterations
    start_index = (i - 1) * window + 1; % Indice di inizio della finestra corrente
    end_index = i * window; % Indice di fine della finestra corrente
    iter = data(start_index:end_index);
    [h(i), p(i), jbstat(i), critval(i)] = jbtest(iter,0.05);
%     if h(i)== 1
%         h(i)="reject"
%     else
%         h(i)="accept"
%     end
end

jb_test = array2table([h p jbstat critval])
rejection_rate = sum(h)/n_iterations
mean(jbstat)
%%
figure
plot(h)
hold on
plot(critval)
histcounts(jbstat,[5.7728 +inf])
%*************************************************************************
%% *ACF ANALYSIS* 
%*************************************************************************
figure('Color',[1 1 1]) %Clusters Plot with Standard Deviation > 2
clusters = data_standardized>2.5 ;
plot(abs(data_standardized.*clusters));%change the value to check for cluster presence
title('Clusters','interpreter','latex')
lags = 1000;
acf_ret = autocorr(data, lags);
acf_abs_ret = autocorr(abs(data), lags);
acf_squared_ret = autocorr(data.^2, lags);

acf_matrix = struct(...
    'Lag', (0:lags)', ...
    'ACF_Returns',acf_ret,...
    'ACF_Absolute_Value_Returns',acf_abs_ret,...
    'ACF_Squared_Returns',acf_squared_ret);
acf_matrix = struct2table(acf_matrix);

%drop first row which is lag0
acf_ret_plot = acf_ret(2:end);
acf_abs_ret_plot = acf_abs_ret(2:end);
acf_squared_ret_plot = acf_squared_ret(2:end);

confidence_interval_pos = ones(size(acf_ret)) *(-1./Descriptive.Number_of_Observations)+2.*(1./sqrt(Descriptive.Number_of_Observations));
confidence_interval_neg = ones(size(acf_ret)) *(-1./Descriptive.Number_of_Observations)-2.*(1./sqrt(Descriptive.Number_of_Observations));

figure('Color',[1 1 1])
hold on
plot(acf_ret_plot, 'blue:', 'LineWidth', 2)
plot(acf_abs_ret_plot, 'black:', 'LineWidth', 2)
plot(acf_squared_ret_plot, 'g:', 'LineWidth',2)
plot(confidence_interval_pos, 'r--', 'LineWidth', 0.4)
plot(confidence_interval_neg, 'r--', 'LineWidth', 0.4)
title('Autocorrelation Functions','interpreter','latex')
legend('Returns', 'Absolute Value', 'Squared','interpreter','latex')
xlabel('Lag','interpreter','latex')
ylabel('Autocorrelation','interpreter','latex')
hold off

%*************************************************************************
%% *BOX AND PIERCE PORTMANTEAU TEST STATISTIC* 
%*************************************************************************
Q_ret = (sum(acf_ret(2:end,:).^2)).*Descriptive.Number_of_Observations;%drop lag 0 again so annoying lol
Q_abs_ret = (sum(acf_abs_ret(2:end,:).^2)).*Descriptive.Number_of_Observations;
Q_square_ret = (sum(acf_squared_ret(2:end,:).^2)).*Descriptive.Number_of_Observations;

CV_Box_Pierce = chi2inv(0.95,lags);

if Q_ret > CV_Box_Pierce
    disp('Returns are autocorrelated')
else
    disp('')
end

if Q_abs_ret > CV_Box_Pierce
    disp('Absolute Returns are autocorrelated')
else
    disp('')
end

if Q_square_ret > CV_Box_Pierce
    disp('Squared Returns are autocorrelated')
else
    disp('')
end

Port_tes = struct(...
    'Portmanteau_Test','Test_Statistics',...
    'Returns',Q_ret,...
    'Absolute_Value_Returns',Q_abs_ret,...
    'Squared_Returns', Q_square_ret,...
    'Critical_Value',CV_Box_Pierce);
Port_tes = struct2table(Port_tes)

%*************************************************************************
%% *LJUNG-BOX TEST STATISTIC (MATLAB VERSION) (check with excel)* 
%*************************************************************************

[h1, pValue1, stat1, cValue1] = lbqtest(data,'lags',lags);
[h2, pValue2, stat2, cValue2] = lbqtest(data.^2,'lags',lags);
[h3, pValue3, stat3, cValue3] = lbqtest(abs(data),'lags',lags);

LJ_Matrix = [h1, pValue1, stat1, cValue1;
               h2, pValue2, stat2, cValue2;
               h3, pValue3, stat3, cValue3];

row_names = {'Returns', 'Squared Returns', 'Absolute Returns'};
column_names = {'hypothesis', 'p-Value', 'Test Statistics', 'Critical Value'};

LJ_Matrix = array2table(LJ_Matrix, 'RowNames', row_names, 'VariableNames', column_names)

%Are Ljung Box Test and Box and Pierce Test equal for large n?
struct('Ljung_Box', stat1,...
   'Box_and_Pierce', Q_ret )

%*************************************************************************
%% *KOLMOGOROV-SMINROFF TEST (for smaller sample sizes)* 
%*************************************************************************
Theoretical_quantiles = sort(data,'ascend');
Empirical_CDF = [1:Descriptive.Number_of_Observations]'./Descriptive.Number_of_Observations;
Theoretical_CDF = normcdf(Theoretical_quantiles,Descriptive.Mean,Descriptive.Standard_Deviation);
difference= abs(Empirical_CDF-Theoretical_CDF);
KOLMOGOROV = max(abs(Empirical_CDF-Theoretical_CDF));

figure('Color',[1 1 1])

subplot(2,2,[1,2]);
hold on;
plot(Theoretical_quantiles, Empirical_CDF);
plot(Theoretical_quantiles, Theoretical_CDF);
xlabel('Theoretical Quantiles','interpreter','latex');
ylabel('CDF Values','interpreter','latex');
title('Empirical vs Theoretical CDF','interpreter','latex');
legend('Empirical CDF', 'Theoretical CDF','interpreter','latex');
hold off;

subplot(2,2,3);
plot(Theoretical_quantiles, difference);
xlabel('Theoretical Quantiles','interpreter','latex');
ylabel('Differences Values','interpreter','latex');
title('Empirical - Theoretical','interpreter','latex');
legend('Difference','interpreter','latex');
grid on;

set(gcf, 'Units', 'normalized', 'Position', [0,0,0.75,0.75]); 
subplot(2,2,[1,2]);
pbaspect([2 1 1]); 

%TO BE COMPLETED
 
%hypothesized_CDF = [Empirical_CDF, Theoretical_CDF];
%[h, p, ksstat, cv] = kstest(data, 'Alpha', 0.05, 'CDF', hypothesized_CDF);

%*************************************************************************
%% *Theoretical vs Actual Probability* 
%*************************************************************************
Theoretical_Prob = [0:0.01:1];

z_Up = norminv((1+Theoretical_Prob)/2);
z_down = -z_Up;
%Now I want to count how many returns are in the interval Z-UP, Z-Down,
%so first you need to standardize the returns, And we extract the returns
%inside the interval, greater than the Lower Bound and Smaller than the
%Upper Bound
OBS_Outside = sum(data_standardized < z_down)+sum(data_standardized > z_Up);
%create dataset for those observations not falling in the interval
%chosen = zeros(Descriptive.Number_of_Observations,1);
%for i = 1:length(data_standardized)
    %if data_standardized(i) > z_Up
        %chosen(i) = data_standardized(i);
    %elseif data_standardized(i) < z_down
        %chosen(i) = data_standardized(i);
    %end
%end
Observed_frequency_Inside = (Descriptive.Number_of_Observations-OBS_Outside)/Descriptive.Number_of_Observations;
Observed_frequency_Outside = OBS_Outside/Descriptive.Number_of_Observations;
Theoretical_Frequency = 1-Theoretical_Prob;

Empirical_vs_Theoretical = struct(...
    'Z_UP', z_Up',...
    'Z_DOWN', z_down',...
    'Observations_Outside',OBS_Outside',...
    'Observed_frequency_Inside',Observed_frequency_Inside',...
    'Theoretical_Probability', Theoretical_Prob');
Empirical_vs_Theoretical=struct2table(Empirical_vs_Theoretical)
%%
%data_standardized = (data-Descriptive_table.Mean)/Descriptive_table.Standard_Deviation;
data_standardized = data;
gg = [Dates(2:end,1) array2table(data_standardized)];
a = Descriptive_table.Mean-3*Descriptive_table.Standard_Deviation;
observed_frequency = ((histcounts(data_standardized,[a -a]))/2515)*100
figure
scatter(gg,"Date","data_standardized","filled",'SizeData',3)
hold on
line(xlim,[a, a], 'Color',"r", 'LineWidth', 1);
line(xlim,[-a, -a], 'Color',"r", 'LineWidth', 1);
title('Portfolio Returns, 0.997 confidence bands','Interpreter','latex')
xlabel('Dates','Interpreter','latex')
ylabel(" ")
hold off
less_than=100*((sum(data<a)+sum(data>-a))/Descriptive_table.Number_of_Observations)
c = ((norminv(1-0.997))*Descriptive_table.Standard_Deviation+Descriptive_table.Mean)