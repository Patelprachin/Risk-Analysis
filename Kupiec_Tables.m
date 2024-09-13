benchmark_PV = varbacktest(data(window:end),-Var_Parametric_in_returns(:,[91,100]),'VaRLevel',[0.9 0.99]);
benchmark_NPV = varbacktest(data(window:end),Var_Non_Parametric_in_returns(:,[91,100]),'VaRLevel',[0.9 0.99]);
benchmark_MCV = varbacktest(data(window:end),-Var_Monte_carlo(:,[91,100]),'VaRLevel',[0.9 0.99]);

violations_PV = cc(benchmark_PV,"TestLevel",0.95);
violations_NPV = cc(benchmark_NPV,"TestLevel",0.95);
violations_MCV = cc(benchmark_MCV,"TestLevel",0.95);

a = [violations_PV.Failures violations_NPV.Failures violations_MCV.Failures]
violations_PV = runtests(benchmark_PV,"TestLevel",0.95);
% violations_NPV = runtests(benchmark_NPV,"TestLevel",0.95);
% violations_MCV = runtests(benchmark_MCV,"TestLevel",0.95);
% [violations_PV.TL violations_NPV.TL violations_MCV.TL]


