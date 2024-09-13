%ZoomInd   = 1:1:(size(data,1))'; uncomment to see the whole big fatass
%plot
ZoomInd   = 120:1:300;
VaRData   = [Var_Parametric_in_returns_plot Var_Non_Parametric_in_returns_plot];
VaRFormat = {'-','--'};
D = ZoomInd;
R = data(ZoomInd);
N = Var_Parametric_in_returns_plot(ZoomInd);
H = Var_Non_Parametric_in_returns_plot(ZoomInd);
IndN95    = (R < -N);
IndHS95   = (R < -H);
figure
bar(D,R,0.5,'FaceColor',[0.7 0.7 0.7]);
hold on
for i = 1 : size(VaRData,2)
    stairs(-VaRData(:,i),VaRFormat{i});
end

ylabel('VaR')
xlabel('Date')
legend({'Returns','VaR Parametric','VaR non Parametric'},'Location','Best','AutoUpdate','Off')
title('95% VaR violations for different models')
ax = gca;
ax.ColorOrderIndex = 1;
plot(D(IndN95),-N(IndN95),'o',D(IndHS95),-H(IndHS95),'o','MarkerSize',5,'LineWidth',1.5)
xlim([D(1)-1, D(end)+1])
hold off;
