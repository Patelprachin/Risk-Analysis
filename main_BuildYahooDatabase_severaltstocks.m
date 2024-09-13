clear all 
close all 
Date = datestr(today());
 
filenameTICKER = 'Tickers_2023.xlsx'
filenameOUTPUT = ['RA_202324_' num2str(year(Date)) '.xlsx']

startdate = '1-Jan-2014'; %starting date
enddate = Date;%Date; %ending date
enddate = '31-Dec-2023';%Date; %ending date

interval = '1d';%sampling frequency

[void tickers] = xlsread(filenameTICKER, 'Tickers', 'A2:A26') 
%tickers = [{'INTC'}, {'JPM'}, {'AA'},{'PG'},{'MSFT'}]'

numtickers = size(tickers); 

for idx = [1:numtickers]% numtickers
    disp('Request Historical Prices')
    symbol = tickers{idx}

    %Download historical prices
    StockData = get_MarketDataViaYahoo(symbol, startdate, enddate, interval);

    %Extract dates and Adj Close
    dates = datenum(StockData.Date);
    adjclose = StockData.AdjClose;
    %Write to Excel file the downloaded data
    writetable(StockData, filenameOUTPUT, 'Sheet', symbol)

    %Compute Log-returns and build a Table
    LogRet = log(adjclose(2:end)./adjclose(1:end-1));
    TabLogret = table(datestr(dates(2:end)),LogRet);
    TabLogret.Properties.VariableNames = {'Dates', 'Log-Ret'};

    %Write to Excel filenameOUTPUT the log-returns
    writetable(TabLogret, filenameOUTPUT, 'Sheet', [symbol '_logret'])

    %Merge data set
    if idx == 1
        merge_date_closing = dates;
        merge_close = adjclose;
        merge_date_logret = dates(2:end);
        merge_logret = LogRet;
        Names{1} = 'Dates';
        Names{2} = symbol;
    end

    if idx>1

        %merge adjusted closing prices
        [merge_date_closing index1 index2] = intersect(merge_date_closing, dates);
        merge_close = [merge_close(index1,:), adjclose(index2)];
        %merge log-returns
        [merge_date_logret index1 index2] = intersect(merge_date_logret, dates(2:end));
        merge_logret = [merge_logret(index1,:), LogRet(index2)];
        Names{idx+1} = symbol;

    end

end

%Write to Excel and store the database
T_closing = splitvars(table(merge_date_closing, merge_close));
T_closing.Properties.VariableNames = Names;
T_closing.Dates=datestr(T_closing.Dates)
writetable(T_closing , filenameOUTPUT, 'Sheet', 'Closing Prices')    

T_logret = splitvars(table(merge_date_logret, merge_logret));
T_logret.Properties.VariableNames = Names;
T_logret.Dates=datestr(T_logret.Dates)
writetable(T_logret, filenameOUTPUT, 'Sheet', 'Log-returns')    

writecell(tickers, filenameOUTPUT, 'Sheet', ['Tickers'])    

