
Varlevel = 100;
Vars = [Var_Parametric_in_returns(:,Varlevel) -Var_Non_Parametric_in_returns(:,Varlevel)...
    Var_Monte_carlo(:,Varlevel)];
clusters= ones(size(Var_Monte_carlo,1),1)*999;
for i = 1:size(Var_Monte_carlo,1)
    for j = 1:size(Vars,2)
        if data(i)<Vars(i,j);
        clusters(i,j)= (data(i)-Vars(i,j))/std(data);
    else
        clusters(i,j)=0;
        end
    end
end

figure('Color',[1 1 1])
subplot(1,3,1)
plot(clusters(:,1),Color=[0 0.4470 0.7410])
xlabel('Dates','Interpreter','latex')
ylabel('Losses','Interpreter','latex')
title('PV volatility clusters','Interpreter','latex')
subplot(1,3,2)
plot(clusters(:,2),Color=[0.8500 0.3250 0.0980])
xlabel('Dates','Interpreter','latex')
ylabel('Losses','Interpreter','latex')
title('NPV volatility clusters','Interpreter','latex')
subplot(1,3,3)
plot(clusters(:,3),Color=[0.9290 0.6940 0.1250])
xlabel('Dates','Interpreter','latex')
ylabel('Losses','Interpreter','latex')
title('MCV volatility clusters','Interpreter','latex')