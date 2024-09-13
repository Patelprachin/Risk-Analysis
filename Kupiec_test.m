N = size(violations_Parametric,1);
x = sum(violations_Parametric)';
p = 1- confidence_interval';
LPOF = -2*((N-x).*log(N.*(1-p)./(N-x))+ x.*log(N*p./x));
p_value = 1-chi2cdf(LPOF,0.1);
significance = chi2cdf(confidence_interval,1)';
for i = 1:size(confidence_interval,2)
        if p_value(i)>=significance(i)
            rejection1(i) = 'reject';
        else
            rejection1(i)= 'accept';
        end
end
benchmark = varbacktest(data(120:end),Var_Parametric_in_returns,'VaRLevel',[linspace(0.01,0.99,100)]);
VaR_tests = runtests(benchmark,"TestLevel",0.95);
POF = pof(benchmark,"TestLevel",0.95);
[POF(:,4) array2table(rejection1')]
a = -2*((N-138).*log(N.*(1-0.9504)./(N-138))+ 138.*log(N*(1-0.9504)./x));