f = @(x) (5 - x)*exp(x) - 5;

a = 4;
b = 5;

tol = 10e-6;

format long;

[root, count] = bisection(a, b, f, tol)

function [x, count] = bisection(a, b, f, tol)
    ntol=ceil(log((b-a)/tol)/log(2));
    for j = 1:ntol
        x = (a+b)/2;
        fx = f(x);
        fa = f(a);
    
        if (fx*fa < 0)
            b = x;
        else
            a = x;
        end
        interval = [a, b]
    end
    count = ntol;
end