%% VaRs Statistical Analysis Plots*****************************************
figure('color', [1 1 1])
plot(confidence_interval, Var_Parametric_in_returns(1,:), 'LineWidth', 1,'LineStyle','-.')
hold on
plot(confidence_interval, -Var_Monte_carlo(1,:), 'LineWidth', 1,'LineStyle','--')
plot(confidence_interval, Var_Non_Parametric_in_returns(1,:), 'LineWidth', 1)
title('VaRs vs Confidence Levels','interpreter','latex')
xlabel('Confidence Level','interpreter','latex')
ylabel('VaRs Value','interpreter','latex')
legend('Parametric Var', 'Monte Carlo VaR','Non-Parametric VaR','Interpreter','latex','Location','Best')

