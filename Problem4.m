n = 1000;

T = 2*pi;
alpha = 0.7*pi;
beta = 0.7*pi;
%alpha = 0.01*pi;
%beta = -0.85*pi;

h = T/n;
t = ((1:n-1)*h).';
tol = 10e-3;

theta0 = alpha*cos(t) + 0.5*sin(t);
%theta0 = 0.7*ones(n-1, 1);
%theta0 = (((beta - alpha) / (2*pi))*t + alpha).*cos(2*t);

tvals = [0 t.' T];
theta0Vals = [alpha theta0.' beta];

plot(tvals, theta0Vals, 'b', "LineWidth", 2)
title("Newton's Method Solution to the Pendulum BVP")
ylim([-5, 5])
hold on;

theta = nDnewton(n, alpha, beta, h, t, theta0, tol);

thetaVals = [alpha theta.' beta];
plot(tvals, thetaVals, 'g', "LineWidth", 2)
legend("\theta^{(0)}", "\theta^{(1)}", "\theta^{(2)}", "\theta^{(3)}", "\theta^{(4)}", "\theta^{(5)}", "\theta^{(6)}")

function theta = nDnewton(n, alpha, beta, h, t, theta0, tol)
    prevTheta = theta0;
    theta = theta0;
    error = norm(prevTheta);
    count = 0;
    tvals = [0 t.' 2*pi];
    while (error > tol)
        J = jacobian(theta, h, n);
        ftheta = f(theta, alpha, beta, h);
        Jinvftheta = J \ ftheta;

        theta = prevTheta - Jinvftheta;

        error = norm(theta - prevTheta);

        prevTheta = theta;
        count = count + 1;
        
        if (error > tol)
            thetaVals = [alpha theta.' beta];
            plot(tvals, thetaVals, "LineWidth", 2)
        end
    end
    count
end

function J = jacobian(theta, h, n)
    mainDiag = (-2/h^2)+cos(theta);
    otherDiags = (1/h^2)*ones(n-2, 1);
    J = diag(mainDiag, 0) + diag(otherDiags, -1) + diag(otherDiags, 1);
end

function ftheta = f(theta, alpha, beta, h)
    ftheta = zeros(length(theta), 1);
    for j = 1:length(theta)
        ftheta(j) = fj(theta, j, alpha, beta, h);
    end
end

function ftheta = fj(theta, j, alpha, beta, h)
    fullTheta = [alpha; theta; beta];
    ftheta = (1/h^2) * (fullTheta(j)-2*fullTheta(j+1)+fullTheta(j+2))+sin(fullTheta(j+1));
end