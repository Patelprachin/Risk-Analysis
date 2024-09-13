clear all
close all
clc 
format short
%***********************************************************************
%% QUESTION 2
%***********************************************************************
Stock_names = ["INTC","JPM","AA","PG","MSFT"];
start_date = ['01-Jan-2014'] ;
end_date = ['31-Dec-2023'] ;
Dates = get_MarketDataViaYahoo('INTC',start_date, end_date,'1d');


df = [Dates(:,1) ];
for i = Stock_names
    df_values = get_MarketDataViaYahoo(convertStringsToChars(i),start_date, end_date,'1d');
    df_adj = df_values(:,6); %Index for col = 6 because it's the 6TH scraped column from YFinance
    df_adj.Properties.VariableNames = [i];
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
ret(1,:)= [];
ret_1 = ret(1:1258,:);
ret_2 = ret(1258:end,:);
%***************************************************************************
%% Determine the number of assets and the sample size
%***************************************************************************
NumTickers = size(Stock_names);
[NumObs NumAssets]= size(ret_1);


%Assign Confidence level
alpha = 0.95;

%Assign weigths
w = ones(NumAssets,1)/NumAssets;

%compute portfolio return
LogRetp = ret_1*w;

%Estimate the Covariance Matrix using the Sample Covariance Matrix
Sigma = cov(ret_1);
%Compute portfolio variance
sg2p = w'*Sigma*w;
%Compute portfolio VaR (assume zero mean)
z = norminv(1-alpha,0,1);
VaR_g = - z* sg2p^0.5;

%Compute Marginal VaR
MVaR_g = - z*Sigma*w/sg2p^0.5;

%Compute Component VaR
CVaR_g = w.*MVaR_g;
chk=[sum(CVaR_g) VaR_g]

%Compute Component VaR in percentage
CVaR_g_p = CVaR_g/sum(CVaR_g);

T_Component_g = table(w, MVaR_g,CVaR_g, CVaR_g_p);
T_Component_g.Properties.VariableNames = {'Weigths' 'MVaR' 'CVaR' 'CVaR\%'};
T_Component_g.Properties.RowNames = Stock_names

%Build the minimum variance portfolio
x0 = ones(NumAssets,1)/NumAssets;
w_mv = fmincon(@(x)  x'*Sigma*x, x0, [], [], ones(1,NumAssets), 1, zeros(NumAssets,1),ones(NumAssets,1)) ;
sg2mv = w_mv'*Sigma*w_mv;
MVaR_mv = - z*Sigma*w_mv/sg2mv^0.5;
CVaR_mv = w_mv.*MVaR_mv;
CVaR_mv_p = CVaR_mv/sum(CVaR_mv);

%Build the risk-parity portfolio
x0 = ones(NumAssets,1)/NumAssets;
w_rp = fmincon(@(x)  std(x.*Sigma*x/(x'*Sigma*x)^0.5), x0, [], [], ones(1,NumAssets), 1, zeros(NumAssets,1),ones(NumAssets,1)) ;
sg2rp = w_rp'*Sigma*w_rp;
MVaR_rp = - z*Sigma*w_rp/sg2rp^0.5;
CVaR_rp = w_rp.*MVaR_rp;
CVaR_rp_p = CVaR_rp/sum(CVaR_rp);

T_Component_rp = table(w_rp, MVaR_rp,CVaR_rp, CVaR_rp_p);
T_Component_rp.Properties.VariableNames = {'Weigths' 'MVaR' 'CVaR' 'CVaR\%'};
T_Component_rp.Properties.RowNames = Stock_names

T_Component_g = table(CVaR_g_p, CVaR_mv_p , CVaR_rp_p );
T_Properties.VariableNames = {"INTC","JPM","AA","PG","MSFT"};
T_Component_g.Properties.RowNames = Stock_names

figure('Color',[1 1 1])
colors = jet(numel(Stock_names));
% Plotting the pie chart for CVaR_g_p
subplot(1,3,1);
pie(CVaR_g_p, Stock_names);
colormap(colors);
title('CVaR general portfolio','Interpreter','latex');

% Plotting the pie chart for CVaR_mv_p
subplot(1,3,2);
pie(CVaR_mv_p, Stock_names);
colormap(colors);
title('CVaR minimum variance portfolio','Interpreter','latex');

% Plotting the pie chart for CVaR_rp_p
subplot(1,3,3);
pie(CVaR_rp_p, Stock_names);
colormap(colors);
title('CVaR risk parity portfolio','Interpreter','latex');
%***********************************************************************
%% risk parity portfolio weights on second half -- returns
%***********************************************************************
RP_port_ret2 = ret_2*w_rp;
EW_port_ret2 =ret_2 * w;
Port_value_EW = prices(1258:end,:)*w;
Port_value_RP = prices(1258:end,:)*w_rp;
%***********************************************************************
%% VaRs
%***********************************************************************
n=1;
confidence_interval = 0.95;
window = 300;
for i = 1:size(RP_port_ret2,1)-window
    data = RP_port_ret2(i:i+window-1);
    mu = mean(data)*n;
    SD = std(data) * sqrt(n);
    z = norminv(1 - confidence_interval);
    VaR_Parametric_loop(i) = -(mu + SD * z);
end
mu = mean(RP_port_ret2) * n; % Sample Mean
SD = std(RP_port_ret2) * sqrt(n);
z = norminv(1 - confidence_interval);
VaR_Parametric = -(mu + SD * z);
VaR_non_Parametric = -prctile(RP_port_ret2,(1-confidence_interval)*100); 
%Number of Violations
for i = 1:size(RP_port_ret2,1)
    if RP_port_ret2(i)> VaR_Parametric
        Parametric_violations(i)=1;
    else
        Parametric_violations(i)=0;
        if RP_port_ret2(i)> VaR_non_Parametric
            Non_Parametric_violations(i)=1;
        else 
            Non_Parametric_violations(i)=0;
        end
    end
end
for i = 1:size(RP_port_ret2,1)
    if RP_port_ret2(i)> VaR_Parametric_loop
        Parametric_violations_loop(i)=1;
    else
        Parametric_violations_loop(i)=0;
    end
end
Measures_1 = struct2table(struct(...
    'Sharpe_Equally_Weighted', sharpe(EW_port_ret2, 4.315/100),...
    'Sharpe_Risk_Parity', sharpe(RP_port_ret2, 4.315/100)));

Measures_2 = struct2table(struct(...
    'Max_Drawdown_Equally_Weighted', maxdrawdown(Port_value_EW),...
    'Max_Drawdown_Risk_Parity', maxdrawdown(Port_value_RP)));


Measures_3 = struct2table(struct(...
    'Parametric_VaR_RP', VaR_Parametric*100,...
    'violations_Parametric', (sum(Parametric_violations)),...
    'Non_Parametric_Var_RP', VaR_non_Parametric*100,...
    'violations_Non_Parametric', (sum(Non_Parametric_violations)),...
    'Parametric_VaR_violations_loop', sum(Parametric_violations_loop)));
%***********************************************************************
%% Plots
%***********************************************************************
figure('Color',[1 1 1])
plot(cumsum(RP_port_ret2))
hold on
plot(cumsum(EW_port_ret2))
% plot(repmat(VaR_Parametric,size(EW_port_ret2,1)))
% plot(repmat(VaR_non_Parametric,size(EW_port_ret2,1)))
legend('Risk Parity Portfolio','Equally Weighted Portfolio', 'Location','Best','AutoUpdate','Off', 'Interpreter','latex')
ylabel('Cumulative Returns','Interpreter','latex')
xlabel('Dates', 'Interpreter','latex')
grid minor
hold off

%*************************************************************************
%% Violations Plot
%*************************************************************************
start_plot = 1;
end_plot = size(VaR_Parametric_loop,2);
ZoomInd= start_plot:1:end_plot; 
%VaRData   = [Var_Parametric_in_returns_plot Var_Non_Parametric_in_returns_plot std_data_plot];
a = repmat(VaR_Parametric,size(EW_port_ret2,1)-window,1);
b = repmat(VaR_non_Parametric,size(EW_port_ret2,1)-window,1);
VaRData   = [VaR_Parametric_loop' a];
Var_Parametric_in_returns_loop_plot = VaR_Parametric_loop;
Var_Parametric_in_returns_plot = a;
VaR_non_Parametric_plot = b;
VaRFormat = {'-','--','-.'};
D = ZoomInd;
R = RP_port_ret2(ZoomInd);
N = Var_Parametric_in_returns_loop_plot(ZoomInd)';
H = Var_Parametric_in_returns_plot(ZoomInd);
I = VaR_non_Parametric_plot(ZoomInd);

IndN95    = (R < -N);
IndHS95   = (R < -H);
INdHZ95  = (R<-I);
figure('color', [1 1 1]);
bar(D,R,0.5,'FaceColor',[0.7 0.7 0.7]);
hold on
plot(-Var_Parametric_in_returns_loop_plot,'Color','b')
plot(-Var_Parametric_in_returns_plot,'Color','r')
plot(-VaR_non_Parametric_plot,'Color','g')
ylabel('VaR','Interpreter','latex')
xlabel('Date','Interpreter','latex')
legend({'Returns','VaR Parametric','VaR non Parametric'},'Location','Best','AutoUpdate','Off', 'Interpreter','latex')
title('95% VaR violations for different models','Interpreter','tex')
ax = gca;
ax.ColorOrderIndex = 1;
plot(D(IndN95),-N(IndN95),'o',D(IndHS95),-H(IndHS95),'o',D(INdHZ95),-I(INdHZ95),'o','MarkerSize',5,'LineWidth',1.5)
xlim([D(1)-1, D(end)+1])
grid minor
hold off;

figure('color',[1 1 1])
numBins = round(sqrt(size(RP_port_ret2,1)));
histEdges = linspace(min(RP_port_ret2), max(RP_port_ret2), numBins);
%histogram(RP_port_ret2, histEdges)
hold on 
q = quantile(RP_port_ret2, 0.05);
q1 = -Var_Parametric_in_returns_plot(1,1);
max_Var_Parametric = max(-Var_Parametric_in_returns_plot(1,1));

returns_left_of_VaR_non_Parametric = RP_port_ret2(RP_port_ret2 < q);
histogram(returns_left_of_VaR_non_Parametric, histEdges,'FaceColor',[0.6 0 0])

returns_between_Var_Parametric_and_max = RP_port_ret2(RP_port_ret2 >= q & RP_port_ret2 <= max_Var_Parametric);
histogram(returns_between_Var_Parametric_and_max, histEdges,'FaceColor',[0.6 0.6 0.6])

returns_right_of_Var_Parametric = RP_port_ret2(RP_port_ret2 > max_Var_Parametric);
histogram(returns_right_of_Var_Parametric, histEdges,'FaceColor',[0.2 0.8 0.2])

line([q, q], ylim, 'Color', [1 0.5 0.3], 'LineWidth', 1, 'Linestyle',':');
line([q1 q1], ylim, 'Color', [1 0.3 0.3], 'LineWidth', 1,'Linestyle','--');

xlabel('Risk Parity Portfolio Returns', 'Interpreter','latex')
ylabel('Frequencies', 'Interpreter','latex')
legend('Risky Zone','Tolerance Zone','Acceptance zone','VaR Non Parametric','VaR Parametric','Interpreter','latex')
grid minor
hold off
