window = 120;
n = round(length(data)/window);

final_value = zeros(n + 1, size(data, 2));
final_value(1, :) = data(1, :);

for t = 1:n
    index = t * window;
    if index <= length(data)
        final_value(t + 1, :) = data(index, :);
    else
        final_value(t + 1, :) = data(end, :); 
    end
end

[final_value([1,2,3],:) data([1,60,120],:)]


