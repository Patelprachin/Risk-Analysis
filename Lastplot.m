rng(123)
y = trnd(5,2515,1);
y1 = randn(2515,1);
%%
figure
subplot(1,2,1)
qqplot(y,data)
subplot(1,2,2)
qqplot(y1,data)
%%
%a = normcdf(data,Descriptive.Mean,Descriptive.Standard_Deviation);
a = normcdf(VaR,mean(data),std(data))
% rng default   % For reproducibility
% x = data;
% mu = Descriptive_table.Mean;
% sigma = Descriptive_table.Standard_Deviation;
% xbar = mean(x);
% s = std(x);
% t = (x-mu)./(sigma/sqrt(size(data,1)));
% % b = tcdf(t,5);
% [a1,a2] = ttest(data,Descriptive.Mean,0.05,'right');
 figure
% histogram(b)
histogram(a)
xlabel('Probability of negative returns','Interpreter','latex')
ylabel('Frequency','Interpreter','latex')
title('Transform Probability','Interpreter','latex')
%%
b = norminv(a, 0, 1);
inf_rows_indices = find(any(isinf(b), 2));
figure('Color',[1 1 1])
histogram(b)
numBins = round(sqrt(size(b,1)));
histEdges = linspace(min(b), max(b), numBins);
%xlabel('Probability of negative returns','Interpreter','latex')
ylabel('Frequency','Interpreter','latex')
title('Normal Random Variable','Interpreter','latex')
%%
X = [ones(size(data,1),1) data];
X(inf_rows_indices,:)= [];
b(inf_rows_indices,:) = [];
betas = regress(b,X)