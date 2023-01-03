%=========================================%
% Boundary value problem (BVP) for        %
% ODE of second order:                    %
% u'' = u' + 2*y + cos(x)                 %
% (1) u(0) = -0.3                         %
% (2) u(pi/2) = -0.1                      %
%=========================================%
% Exact solution for comparison:          %
% u(x) = -1/10 * (sin(x) + 3*cos(x))      %
%=========================================%
% Numerical method: conservative finite   %
% difference method.                      %
%=========================================%

clear; clc;
tic;

%=================%
% Spatial grid.   %
%=================%
a = 0; b = pi/2; 
x = linspace(a,b,200);
h = x(2) - x(1);
N = length(x);

%===================================%
% Tabulating the exact solution.    %
%===================================%
ur = zeros(1,N);
for i = 1 : N
    ur(i) = -1/10 * (sin(x(i)) + ...
        3 * cos(x(i)));
end

%========================%
% Approximate solution.  %
%========================%
y = zeros(1,N);

%==============================%
% Left boundary condition.     %
%==============================%
y(1) = ur(1);

%=============================%
% Right boundary condition.   %
%=============================%
y(N) = ur(N);

%=================================%
% Matrix of the discrete linear   %
% system of equations.            %
%=================================%
system_matrix = zeros(N,N);

%============================%
% Diagonals of the system.   %
%============================%
% Lower diagonal.            %
%============================%
adiag = zeros(1,N);

%============================%
% Main diagonal.             %
%============================%
bdiag = zeros(1,N);

%============================%
% Upper diagonal.            %
%============================%
cdiag = zeros(1,N);

%============================%
% Right side of the system.  %
%============================%
right_side = zeros(N,1);

%===================================%
% Filling with values adiag.        %
%===================================%
adiag(1) = 0;
for ii = 2 : N-1
    adiag(ii) = 1 + h/2;
end
adiag(end) = 0;

%===================================%
% Filling with values bdiag.        %
%===================================%
bdiag(1) = 1;
for ii = 2 : N-1
    bdiag(ii) = -2*(1+h^2);
end
bdiag(end) = 1;

%===================================%
% Filling with values cdiag.        %
%===================================%
cdiag(1) = 0;
for ii = 2 : N-1
    cdiag(ii) = 1 - h/2;
end
cdiag(end) = 0;

%======================================%
% Filling with values the right side.  %
%======================================%
right_side(1) = ur(1);
for ii = 2 : N-1
    right_side(ii) = h^2 * cos(x(ii));
end
right_side(end) = ur(end);

%====================================%
% Defining the matrix of the linear  %
% system of equations.               %
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

%===================================%
% Progonka method for solving the   %
% linear system of equations.       %
%===================================%
y = Progon(system_matrix, right_side)';

%===========%
% Plot      %
%===========%
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
       '\bf{of second order with Dirichlet }', ...
       '\bf{boundary condition}'])

%==================%
% Absolute error   %
%==================%
err = abs(y - ur);
display(['Maximal error: ', num2str(max(err))])
time = toc;
display(['Elapsed time: ', num2str(time)])
