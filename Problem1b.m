f = @(x) (5 - x)*exp(x) - 5;
fd = @(x) (4 - x)*exp(x);

x0 = 5;

tol = 10e-8;

[root, count] = newton(x0, f, fd, tol)

function [x, count] = newton(x0, f, fd, tol)
    x = x0;
    fx = f(x);
    count = 0;
    while (abs(fx) > tol)
        fdx = fd(x);
        x = x - (fx / fdx)
        count = count + 1;
        fx = f(x)
    end
end