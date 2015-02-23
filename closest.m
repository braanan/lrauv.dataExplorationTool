function [index] = closest(x,y)

% Find closest value of x in y vector

index = zeros(size(x));
for j = 1:length(x)
    [~,index(j)] = min(abs(abs(x(j))-y));
end
end