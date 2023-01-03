%==========================================================%
% Conservative difference scheme                           %
% for ODE of second order with                             % 
% volatile coefficients and                                %
% Dirichlet boundary conditions:                           %
% (k(x)u'(x))' - q(x)u(x) = -f(x), 0 < x < 1,              %
% u(0) = 0                                                 %
% u(1) = sin(1) - 1 = -0.1585                              %
%==========================================================%
% Exact solution for comparison:                           %
% u(x) = sin(x) - x^2                                      %
%==========================================================%
% Data:                                                    %
% k(x) = x^2+1                                             %
% q(x) = x                                                 %
% f(x) = -x^3 + 6*x^2 + (x^2+x+1)*sin(x) - 2*x*cos(x) + 2  %
%==========================================================%

clear; clc; tic;

%===========%
% Grid.     %
%===========%
a = 0; b = 1;
x = linspace(a,b,200);
h = x(2) - x(1);
N = length(x);

%========================%
% Approximate solution.  %
%========================%
y = zeros(1,N);

%===================================%
% Tabulating the exact solution.    %
%===================================%
ur = zeros(1,N);
ur = sin(x) - x.^2;

%===========================%
% Left boundary condition.  %
%===========================%
y(1) = ur(1);

%===========================%
% Right boundary condition. %
%===========================%
y(N) = ur(N);

%===================%
% Function k(x).    %
%===================%
kF = @(x) x.^2+1;

%==================%
% Function q(x).   %
%==================%
qF = @(x) x;

%==================%
% Function f(x).   %
%==================%
fF = @(x) -x.^3 + 6 * x.^2 + ...
          (x.^2 + x + 1) .* sin(x) - ...
          2 * x .* cos(x) + 2;
          
%===============================%
% Matrix of the linear system.  %
%===============================%
system_matrix = zeros(N,N);

%========================%
% Lower diagonal adiag.  %
%========================%
adiag = zeros(1,N);
adiag(1) = 0;
for ii = 2 : N-1
    adiag(ii) = kF(x(ii) - h/2);
end
adiag(end) = 0;

%========================%
% Main diagonal bdiag.   %
%========================%
bdiag = zeros(1,N);
bdiag(1) = 1;
for ii = 2 : N-1
    bdiag(ii) = -(kF(x(ii) + h/2) + ...
        kF(x(ii) - h/2) + h^2 * qF(x(ii)));
end
bdiag(end) = 1;

%========================%
% Upper diagonal cdiag.  %
%========================%
cdiag = zeros(1,N);
cdiag(1) = 0;
for ii = 2 : N-1
    cdiag(ii) = kF(x(ii)+h/2);
end
cdiag(end) = 0;

%=================================%
% Right-hand side of the linear   %
% system of equations.            %
%=================================%
right_side = zeros(N,1);
right_side(1) = y(1);
for ii = 2 : N-1
    right_side(ii) = -h^2 * fF(x(ii));
end
right_side(end) = y(end);

%=====================================%
% Defining the matrix of the linear   %
% system of difference equations.     %
%=====================================%
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

%================================%
% Solving by Progonka method.    %
%================================%
y = Progon(system_matrix, right_side)';

%===================%
% Absolute error.   %
%===================%
err = abs(y - ur);

%========%
% Plot.  %
%========%
figure(1)
% % %
subplot(1,2,1)
% % %
plot(x,y,'b','LineWidth',3)
hold on
grid on
plot(x,ur,'cyan--','LineWidth',3)
set(gca,'FontSize',14)
xlabel('\bf{x}')
ylabel('\bf{y}')
legend('\it{Approximated solution}', ...
       '\it{Exact solution}')
% % %
subplot(1,2,2)
% % %
plot(x,err,'r:','LineWidth',3)
hold on, grid on
set(gca,'FontSize',14)
xlabel('\bf{x}')
ylabel('\bf{y}')
legend('\it{Absolute error}')

%==================================%
% Printing important information.  %
%==================================%
display(['Maximal error: ', num2str(max(err))])
time = toc;
display(['Elapsed time: ', num2str(time)])
