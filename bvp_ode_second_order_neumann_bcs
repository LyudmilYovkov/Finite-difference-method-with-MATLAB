%============================================%
% BVP for ODE of 2nd order                   %
% u'' - u = x, 0 < x < 1                     %
% u(0) = 0                                   %
% u'(1) = (e^2-e+1) / (e^2)                  %
%============================================%
% Exact solution for comparison:             %
% u(x) = e^(-x) * (-1 + e^(2*x) - x*e^(x)).  %
%============================================%

clear;
clc;

tic;
%=================%
% Spatial grid    %
%=================%
a = 0; b = 1;
x = linspace(a,b,200);
h = x(2) - x(1);
N = length(x);

%=================================%
% Tabulating the exact solution.  %
%=================================%
ur = ...
    exp(-x) .* (-1 + exp(2*x) - x.*exp(x));

%========================%
% Approximate solution.  %
%========================%
y = zeros(1,N);

%===================================
% Initial condition
%===================================
%y(1) = 0;

%================================%
% Matrix of the discrete linear  %
% system of equations.           %
%================================%
system_matrix = zeros(N,N);

%=================================%
% Diagonals of the sparse system  %
% matrix.                         % 
%=================================%
% Lower diagonal adiag.           %
%=================================%
adiag = zeros(1,N);
adiag(1) = 0;
for ii = 2 : N-1
    adiag(ii) = 1;
end
adiag(end) = -1/h;

%=================================%
% Main diagonal bdiag.            %
%=================================%
bdiag = zeros(1,N);
bdiag(1) = 1;
for ii = 2 : N-1
    bdiag(ii) = -(h^2+2);
end
bdiag(end) = 1/h + h/2;

%=================================%
% Upper diagonal cdiag.           % 
%=================================%
cdiag = zeros(1,N);
cdiag(1) = 0;
for ii = 2 : N-1
    cdiag(ii) = 1;
end
cdiag(end) = 0;

%==================================%
% Right-hand side of the system.   %
%==================================%
alpha = (exp(1)^2 - exp(1) + 1) / (exp(1));
right_side = zeros(N,1);
right_side(1) = 0;
for ii = 2 : N-1
    right_side(ii) = h^2 * x(ii);
end
right_side(end) = alpha - h/2 * x(end);

%====================================%
% Defining the matrix of the system. %
%====================================%
for i0 = 1 : N
    for j0 = 1 : N
        if(i0 == j0+1)
            system_matrix(i0,j0) = adiag(i0);
        end
        if(j0 == i0+1)
            system_matrix(i0,j0) = cdiag(i0);
        end
        if(i0 == j0)
            system_matrix(i0,j0) = bdiag(i0);
        end
    end
end

%===================%
% Progonka method.  %
%===================%
y = Progon(system_matrix, right_side)';

%==========%
% Plot.    %
%==========%
figure(1)
plot(x,y,'b','LineWidth',3)
hold on
grid on
plot(x,ur,'g--','LineWidth',3)
xlabel('\bf{x}')
ylabel('\bf{Solution}')
legend('\bf{Approximated solution}', ...
       '\bf{Exact solution}')
title(['\bf{Boundary value problem for ODE }', ...
       '\bf{of second order with a Neumann }', ...
       '\bf{boundary condition}'])
       
%===================%
% Absolute error.   %
%===================%
err = abs(y - ur);
display(['Maximal error: ', num2str(max(err))])
time = toc;
display(['Elapsed time: ', num2str(time)])
