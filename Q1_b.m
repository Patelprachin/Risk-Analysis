close all
clc
%*************************************************************************
%% GAUSSIAN-PARAMETRIC APPROACH VaR VS NON-PARAMETRIC APPROACH
%*************************************************************************
window =120;
n = 1; % Time Horizon
confidence_interval = linspace(0.01,0.99,100);
n_iterations = length(data) - window + 1;
%df = zeros(window, n_iterations, numel(confidence_interval));
Var_Parametric_in_returns = zeros(n_iterations, numel(confidence_interval));
Var_Non_Parametric_in_returns = zeros(n_iterations, numel(confidence_interval));
Var_Monte_carlo = zeros(n_iterations, numel(confidence_interval));
Exp_Shortfall_Param = zeros(1,size(confidence_interval,2));
tic
for j = 1:numel(confidence_interval)
    for i = 1:n_iterations
        iterations = data(i:i+window-1);
        %df(:, i, j) = iterations;
        mu = mean(iterations) * n; % Sample Mean
        SD = std(iterations) * sqrt(n);
        %standard_deviation_data(i) = std(iterations);
        z = norminv(1 - confidence_interval(j));
        VaR = -(mu + SD * z);
        Var_Parametric_in_returns(i, j) = VaR;
        %Exp_Shortfall_Param(j) = -(mu-(SD*normpdf(z)/(1-confidence_interval(j))));
        %1st approach to non parametroc var with matlab interpolation
        iter = sort(iterations,'ascend');
        VaR1 = prctile(iter,(1-confidence_interval(j))*100); 
        Var_Non_Parametric_in_returns(i,j) = -VaR1;
        k = mu*window - (((SD*sqrt(window))^2)/2);
        Var_Monte_carlo(i,j) = prctile(normrnd(0,1,500,1) * sqrt(n/250) * ...
        (SD * sqrt(window)) + k * n/250, (1-confidence_interval(j))*100);
    end 
end
toc
%*************************************************************************
%% Plots
%*************************************************************************
VaRlevel = 90;
additional_rows = zeros(window,1);
Var_Parametric_in_returns_plot= [additional_rows; Var_Parametric_in_returns(:,VaRlevel)];%comment to not have plot
Var_Non_Parametric_in_returns_plot= [additional_rows; Var_Non_Parametric_in_returns(:,VaRlevel)];
Var_Monte_carlo_plot = [additional_rows; Var_Monte_carlo(:,VaRlevel)];
%std_data_plot = [additional_rows ; standard_deviation_data'];
a = -Var_Parametric_in_returns_plot;
b = -Var_Non_Parametric_in_returns_plot;
c = Var_Monte_carlo_plot;
%% VaRs Over TimePlots*****************************************************
figure('color', [1 1 1])
plot(data,LineWidth=0.00025)
hold on 
xlim([0 size(data,1)])
plot(a,LineWidth=1.5,Color = [0.8500 0.3250 0.0980],LineStyle="-.")
plot(b,LineWidth=1.5,Color = [0.9290 0.6940 0.1250],LineStyle="--")
plot(c,LineWidth = 1, Color=[0.4660 0.6740 0.1880],LineStyle=":")
hold off
legend('Portfolio Returns','VaR Parametric', 'Var non Parametric','Var MonteCarlo', 'interpreter','latex')
title('VaRs Over Time', 'Interpreter','latex')
xlabel('Dates', 'Interpreter','latex')
ylabel('VaRs vs Returns', 'Interpreter','latex')
grid minor
%% VaRs Statistical Analysis Plots*****************************************
figure('color', [1 1 1])
subplot(3,3,1)
plot(confidence_interval, Var_Parametric_in_returns(1,:), 'LineWidth', 0.5)
title('Parametric VaR vs Confidence Levels','interpreter','latex')
xlabel('Confidence Level','interpreter','latex')
ylabel('Parametric VaR','interpreter','latex')
grid on

subplot(3,3,2)
plot(confidence_interval, -Var_Monte_carlo(1,:), 'LineWidth', 0.5)
title('Monte Carlo VaR vs Confidence Levels','interpreter','latex')
xlabel('Confidence Level','interpreter','latex')
ylabel('Monte Carlo VaR','interpreter','latex')
grid on


subplot(3,3,3)
plot(confidence_interval, Var_Non_Parametric_in_returns(1,:), 'LineWidth', 0.5)
%hold on 
%plot(confidence_interval, Var_Non_Parametric_in_returns2(1,:), 'LineWidth', 0.5)
title('Non Parametric VaR vs Confidence Levels','interpreter','latex')
xlabel('Confidence Level','interpreter','latex')
ylabel('Non Parametric VaR','interpreter','latex')
grid on

subplot(3,3,4)
yyaxis left
plot(Var_Parametric_in_returns_plot(:,1), 'LineWidth', 0.5)
ylabel('Parametric VaR','interpreter','latex')
yyaxis right
plot(data, 'LineWidth', 0.5)
ylabel('Portfolio Returns','interpreter','latex')
title('VaR Over Time (Parametric)','interpreter','latex')
xlabel('Dates','interpreter','latex')
grid on


subplot(3,3,5)
yyaxis left
plot(Var_Monte_carlo(:,1), 'LineWidth', 0.5)
ylabel('Monte Carlo VaR','interpreter','latex')
yyaxis right
plot(data, 'LineWidth', 0.5)
ylabel('Portfolio Returns','interpreter','latex')
title('VaR Over Time (Monte Carlo)','interpreter','latex')
xlabel('Dates','interpreter','latex')
grid on

subplot(3,3,6)
yyaxis left
plot(Var_Non_Parametric_in_returns_plot(:,1), 'LineWidth', 0.5)
ylabel('Non Parametric VaR','interpreter','latex')
yyaxis right
plot(data, 'LineWidth', 0.5)
ylabel('Portfolio Returns','interpreter','latex')
title('VaR Over Time (Non-Parametric)','interpreter','latex')
xlabel('Dates','interpreter','latex')
grid on

subplot(3,3,[7 8])
plot(Var_Non_Parametric_in_returns_plot(:,1),'Linewidth',0.5)
hold on
xlim([0 size(data,1)])
plot(Var_Parametric_in_returns_plot,'LineWidth',0.5)
plot(-Var_Monte_carlo_plot,'LineWidth',0.5)
title('VaRs Over Time','interpreter','latex')
legend('Var Non Parametric','VaR Parametric','VaR Monte Carlo', 'Interpreter', 'latex')
xlabel('Dates','Interpreter','latex')
ylabel('VaRs','Interpreter','latex')
grid minor
hold off

general_mean = (mean(Var_Monte_carlo(:,VaRlevel)) +...
    mean(Var_Parametric_in_returns(:,VaRlevel)) ...
    +mean(Var_Non_Parametric_in_returns(:,VaRlevel)))/3;
subplot(3,3,9)
scatter3(Var_Non_Parametric_in_returns_plot(:,1), Var_Parametric_in_returns_plot,...
    -Var_Monte_carlo_plot, 'b','.');
ylabel('Parametric VaR','interpreter','latex')
xlabel('Non Parametric VaR','interpreter','latex')
zlabel('Monte Carlo VaR', 'Interpreter', 'latex')
title('Scatter Plot','interpreter','latex')
grid minor
%*************************************************************************
%% Metrics
%*************************************************************************
for i = 1:size(Var_Non_Parametric_in_returns,1)
    for j = 1:size(Var_Non_Parametric_in_returns,2)
        if data(i) <= -Var_Non_Parametric_in_returns(i,j)
            violations_Non_Parametric(i,j) = 1;
        else
            violations_Non_Parametric(i,j) = 0;
        end
    end
end

for i = 1:size(Var_Parametric_in_returns,1)
    for j = 1:size(Var_Parametric_in_returns,2)
        if data(i) <= -Var_Parametric_in_returns(i,j)
            violations_Parametric(i,j) = 1;
        else
            violations_Parametric(i,j) = 0;
        end
    end
end
Hit_Rate = 1-confidence_interval;
%*************************************************************************
%% Kupiec-VaR Parametric
%*************************************************************************
benchmark = varbacktest(data(window:end),Var_Parametric_in_returns,'VaRLevel',[linspace(0.01,0.99,100)]);
VaR_tests_parametric = runtests(benchmark,"TestLevel",0.95);
Christoffersen_parametric = cc(benchmark,"TestLevel",0.95);
POF_parametric = pof(benchmark,"TestLevel",0.99);
%*************************************************************************
%% Kupiec-VaR Non Parametric
%*************************************************************************
benchmark = varbacktest(data(window:end),Var_Non_Parametric_in_returns,'VaRLevel',[linspace(0.01,0.99,100)]);
VaR_tests_non_parametric = runtests(benchmark,"TestLevel",0.95);
Christoffersen_non_parametric = cc(benchmark,"TestLevel",0.95);
POF_non_parametric = pof(benchmark,"TestLevel",0.99);
%*************************************************************************
%% Kupiec-VaR Monte Carlo
%*************************************************************************
benchmark = varbacktest(data(window:end),Var_Monte_carlo,'VaRLevel',[linspace(0.01,0.99,100)]);
VaR_tests_Monte_carlo = runtests(benchmark,"TestLevel",0.95);
Christoffersen_Monte_carlo = cc(benchmark,"TestLevel",0.95);
POF_Monte_carlo = pof(benchmark,"TestLevel",0.99);
Table_kupiec = [VaR_tests_parametric(96,:) ;VaR_tests_non_parametric(96,:); VaR_tests_Monte_carlo(96,:)]
%*************************************************************************
%% Violations Plot
%*************************************************************************
start_plot = 1450;
end_plot = 1850;
ZoomInd= start_plot:1:end_plot; 
%VaRData   = [Var_Parametric_in_returns_plot Var_Non_Parametric_in_returns_plot std_data_plot];
VaRData   = [Var_Parametric_in_returns_plot Var_Non_Parametric_in_returns_plot -Var_Monte_carlo_plot];
VaRFormat = {'-','--','-.','---'};
D = ZoomInd;
R = data(ZoomInd);
N = Var_Parametric_in_returns_plot(ZoomInd);
H = Var_Non_Parametric_in_returns_plot(ZoomInd);
E = -Var_Monte_carlo_plot(ZoomInd);
%S = std_data_plot(ZoomInd);
IndN95    = (R < -N);
IndHS95   = (R < -H);
IndEWMA95 = (R < -E);
figure('color', [1 1 1])
bar(D,R,0.5,'FaceColor',[0.7 0.7 0.7]);
hold on
for i = 1 : size(VaRData,2)
    stairs(-VaRData(:,i),VaRFormat{i});
end

ylabel('VaRs','Interpreter','latex')
xlabel('Date','Interpreter','latex')
legend({'Returns','VaR Parametric','VaR non Parametric','VaR MonteCarlo'},'Location','Best','AutoUpdate','Off', 'Interpreter','latex')
title('90% VaR violations for different models','Interpreter','tex')
ax = gca;
ax.ColorOrderIndex = 1;
plot(D(IndN95),-N(IndN95),'o',D(IndHS95),-H(IndHS95),'o',D(IndEWMA95),-E(IndEWMA95),'o','MarkerSize',5,'LineWidth',1.5)
xlim([D(1)-1, D(end)+1])
hold off;