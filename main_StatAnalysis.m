close all
clear all

%Load log-returns and dates
filename = 'Sp500.xlsx';
symbol = 'SP500';
DataSet = readtable(filename, 'Sheet', [symbol '_logret'], 'Range', 'A1')
Dates = DataSet{:,1};
LogRet = DataSet{:,2};

%Basic Statistics: daily returns
nobs = length(LogRet) %sample size
mean_r = mean(LogRet); %sample mean of daily returns
var_r = var(LogRet) %variance of daily returns 
sd_r = var_r^0.5 %standard deviation 
min_r = min(LogRet) % smallest return 
max_r = max(LogRet) % largest return 
sk_r = skewness(LogRet) % sample skewness
k_r = kurtosis(LogRet) %sample kurtosis

%Create a Table to illustrate the results
Synthesis = table(nobs, mean_r, var_r, sd_r, min_r, max_r,sk_r, k_r )
Synthesis.Properties.VariableNames = {'N. Obs', 'Mean', 'Var', 'St.Dev.', 'Min', 'Max', 'Skew', 'Kurt'}
Synthesis.Properties.RowNames = {symbol}

%Compare Gaussian distr. of Histogram
h = figure('Color',[1 1 1])
subplot(2,2,[1:2])
nbins = round(nobs^0.5);
histogram(LogRet, nbins, 'Normalization', 'pdf')
hold on
fplot(@(x) normpdf(x, mean_r, sd_r), [min_r max_r]) 
xlabel('log-returns','interpreter','latex')
title([symbol ': ' 'Empirical vs Gaussian Density'],'interpreter','latex')
xlim([mean_r-6*sd_r mean_r+6*sd_r])

subplot(2,2,3)
histogram(LogRet, nbins, 'Normalization','pdf')
hold on
fplot(@(x) normpdf(x, mean_r, sd_r), [min_r max_r]) 
xlim([[prctile(LogRet,1) prctile(LogRet,45)]])%left tail
xlabel('log-returns','interpreter','latex')
title([symbol ': ' ' Left Tail'],'interpreter','latex')

subplot(2,2,4)
histogram(LogRet, nbins, 'Normalization','pdf')
hold on
fplot(@(x) normpdf(x, mean_r, sd_r), [prctile(LogRet,60) prctile(LogRet,99)]) 
xlim([[prctile(LogRet,55) prctile(LogRet,99)]])%right tail
xlabel('log-returns','interpreter','latex')
title([symbol ': ' ' Right Tail'],'interpreter','latex')
print(h, [symbol '_empvsgauss'],'-dpng')

%QQPlot
Z=(LogRet-mean_r)/sd_r; %standardize returns
thquantile = icdf('normal',[0.5:nobs-0.5]/nobs,0,1)';%theoretical quantiles
emquantile = sort(Z);%empirical quantiles
h=figure('Color',[1 1 1])
plot(thquantile , [thquantile , emquantile ])
axis equal
grid on
xlabel('Theoretical Gaussian quantiles','interpreter','latex')
ylabel('Sample quantiles','interpreter','latex')
print(h, [symbol '_qqplot'],'-dpng')

%Or use directly the Matlab function
h=figure('Color',[1 1 1]);
qqplot(Z)

%Perform Jarque-Bera test
JB = (nobs/6)*sk_r^2+ (nobs/24)*(k_r-3)^2
CV = icdf('chi', 0.95, 2)
if JB>CV
    disp('Reject null of normality')
else
    disp('Accept null of normality')
end

%Using Matlab function
[h,pvalue,jbstat,critval] = jbtest(Z)

%Do losses cluster?
Extreme = zeros(nobs,1);
Extreme(find(Z>3))  = Z(find(Z>3))
Extreme(find(Z<-3))  = Z(find(Z<-3));
h=figure('Color',[1 1 1])
plot(abs(Extreme))
xlabel('Days','interpreter','latex')
ylabel('Losses larger than 3 st. dev.','interpreter','latex')
title('Clustering of extreme losses/gains','interpreter','latex')
print(h, [symbol '_extreme'],'-dpng')

%Are returns autocorrelated?
maxlags = 20;
[acf, lags, bounds, h] = autocorr(LogRet, maxlags )
h=figure('Color',[1 1 1])
%autocorr(Z, maxlags )
xlim([0.5, maxlags+0.5])
ylim([min(acf(2:end)) max(acf(2:end))])
title('Autocorrelation of log-returns','interpreter','latex')
xlabel('Lag','interpreter','latex')
ylabel('Sample autocorrelation','interpreter','latex')
print(h,[symbol '_acf'],'-dpng')

%Are squared returns autocorrelated?
[acfSq, lags, bounds, h] = autocorr(LogRet.^2, maxlags )
h=figure('Color',[1 1 1])
%autocorr(LogRet.^2, maxlags )
xlim([0.5, maxlags+0.5])
ylim([min(acfSq(2:end)) max(acfSq(2:end))])
xlabel('Lag','interpreter','latex')
ylabel('Sample autocorrelation','interpreter','latex')
title('Autocorrelation of squared log-returns','interpreter','latex')
print(h,[symbol '_acf2'],'-dpng')

%Testing the significance of Autocorrelation
BoxPierce  = sum(acf(2:end).^2)*nobs%
LjungBox = sum(acf(2:end).^2./(nobs-[1:maxlags]'))*nobs*(nobs+2)
[h, pValue, stat, cValue] = lbqtest(LogRet)%matlab

BoxPierce2  = sum(acfSq(2:end).^2)*nobs;
LjungBox2 = sum(acfSq(2:end).^2./(nobs-[1:maxlags]'))*nobs*(nobs+2)
[h,pValue,stat2,cValue] = lbqtest(LogRet.^2);%matlab

CV = icdf('chisquare', 0.95, maxlags)
if LjungBox>CV
    disp('There is evidence of serial correlation in log-returns')
else
    disp('There is no evidence of serial correlation in log-returns')
end

if LjungBox2  > CV
    disp('There is evidence of serial correlation in squared log-returns')
else
    disp('There is no evidence of serial correlation in squared log-returns')
end

output = table(lags,acf,acfSq);
output.Properties.VariableNames = {'Lags', 'ACF Ret', 'ACF Ret^2'}
Synthesis.Properties.RowNames = {symbol};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Assignment: Compare theoretical and empirical frequencies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create standardized returns
Z=(LogRet - mean(LogRet))/std(LogRet);

%assign intervals
Range=[0 0.25;
       0.25 0.5;
       0.5 1;
       1 1.5;
       1.5 2;
       2 3;
       3, 100];
for j=1:7 %loop over rows of range   
    TheorPr(j,1) = 2*(normcdf(Range(j,2))-normcdf(Range(j,1)));
    EmpFr(j,1) = sum((abs(Z)<Range(j,2)).*(abs(Z)>Range(j,1)))/nobs
end
%create table
[Range EmpFr TheorPr] 


Xrange=[1:nobs];
h=figure('Color',[1 1 1])
subplot(1,2,1)
plot(datenum(Dates), Z,'.')
hold on
plot(datenum(Dates(Z>3)), Z(Z>3),'r.')
plot(datenum(Dates(Z<-3)), Z(Z<-3),'r.')
yline(3)
yline(-3)
xlabel('Time', 'interpreter','latex')
ylabel('Standardized returns', 'interpreter','latex')
xlim([datenum(Dates(2)) datenum(Dates(end))])
dateaxis('x', 12, Dates(2)) 

subplot(1,2,2)
X = categorical(Range);
bar(X(:,1), [EmpFr*100 TheorPr*100]) 
ylabel('Frequencies (\%)', 'interpreter','latex')
legend('Empirical Frequency','Theoretical Freq.','interpreter','latex')
print(h,[symbol '_emp_vs_theor'],'-dpng')
