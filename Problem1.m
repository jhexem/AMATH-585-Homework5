f = @(x) (5 - x).*exp(x) - 5;

x = linspace(0, 5, 1000);
fx = f(x);

plot(x, fx, 'r', "LineWidth", 2)
title("Plot of f(x)=(5-x)exp(x) - 5")