window =120;
n = 1; % Time Horizon
%confidence_interval = linspace(0.01, 0.99, numel(data) - window - 1); % Uncomment for Cool Plot
confidence_interval = linspace(0.01,0.99,100);
n_iterations = length(data) - window + 1;
%df = zeros(window, n_iterations, numel(confidence_interval));
Var_Parametric_in_returns = zeros(n_iterations, numel(confidence_interval));
Var_Non_Parametric_in_returns = zeros(n_iterations, numel(confidence_interval)); % Initialize Var_Non_Parametric_in_returns
%Var_Non_Parametric_in_returns2 = zeros(n_iterations, numel(confidence_interval));
Exp_Shortfall_Param = zeros(1,size(confidence_interval,2));
tic
for j = 1:numel(confidence_interval)
    for i = 1:n_iterations
        iterations = data(i:i+window-1);
        %df(:, i, j) = iterations;
        mu = mean(iterations) * n; % Sample Mean
        SD = std(iterations) * n;
        z = norminv(1 - confidence_interval(j));
        VaR = -(mu + SD * z);
        Var_Parametric_in_returns(i, j) = VaR;
        Exp_Shortfall_Param(j) = -(mu-(SD*normpdf(z)/(1-confidence_interval(j))));
        %1st approach to non parametroc var with matlab interpolation
        iter = sort(iterations,'ascend');
        VaR1 = prctile(iter,(1-confidence_interval(j))); 
        Var_Non_Parametric_in_returns(i,j) = -VaR1;
    end
end
toc
%%
additional_rows = zeros(window,1);
Var_Parametric_in_returns_plot= [additional_rows; Var_Parametric_in_returns(:,90)];%comment to not have plot
Var_Non_Parametric_in_returns_plot= [additional_rows; Var_Non_Parametric_in_returns(:,90)];
a = -Var_Parametric_in_returns_plot;
b = -Var_Non_Parametric_in_returns_plot;
figure
plot(a,LineWidth=2)
hold on 
plot(b,LineWidth=2)
plot(data,LineWidth=0.025)
hold off
legend('VaR Parametric', 'Var non Parametric', 'Portfolio Returns','interpreter','latex')
title('VaRs Over Time', 'Interpreter','latex')
xlabel('Dates', 'Interpreter','latex')
ylabel('VaRs vs Returns', 'Interpreter','latex')
grid minor