f = @(x) (5 - x)*exp(x) - 5;

x0 = 4;
x1 = 5;

tol = 10e-8;

[root, count] = secant(x0, x1, f, tol)

function [x, count] = secant(x0, x1, f, tol)
    xprev = x0;
    x = x1;
    fx = f(x);
    count = 0;
    while (abs(fx) > tol)
        fxprev = f(xprev);
        xnext = x - ((x - xprev) / (fx - fxprev))*fx;
        xprev = x;
        x = xnext
        fx = f(x)
        count = count + 1;
    end
end