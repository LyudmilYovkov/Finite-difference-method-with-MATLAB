%========================================================%
% Boundary value problem of third type for the           %
% stationery heat equation with volatile coefficients.   %
%========================================================%
% (k(x)u'(x)) - q(x)u(x) = -f(x), 0 < x < X,             % 
% +k(0)u'(0) = beta_1 * u(0) - mu_1,                     %
% -k(1)u'(1) = beta_2 * u(1) - mu_2,                     %
% beta_1 >= 0, beta_2 >= 0, beta_1 + beta_2 > 0,         %
% k(x) >= const > 0, q(x) >= 0                           %
%========================================================%
% Exact solution for comparison: u(x) = x^4 + x^2 - 1    %
%========================================================%
% Data:                                                  %
% k(x) = x^2 + 1 > 0                                     %
% q(x) = (x - 1)^2 >= 0                                  % 
% f(x) = x^6 - 2*x^5 - 18*x^4 - 2*x^3 - 18*x^2 + 2*x - 3 %
% beta_1 = 1, beta_2 = 1, mu_1 = -1, mu_2 = 13, X = 1    % 
%========================================================%

clear; clc;
tic;

% Grid. %
x0 = 0; X = 1; 
x = linspace(x0,X,200);
h = x(2) - x(1);
N = length(x);

% Constants. %
beta_1 = 1; beta_2 = 1;
mu_1 = -1; mu_2 = 13;

% Coefficients. %
k = @(x) x.^2 + 1;
q = @(x) (x - 1).^2;
f = @(x) x.^6 - 2 * x.^5 - 18 * x.^4 - ...
         2 * x.^3 - 18 * x.^2 + 2 * x - 3;

% Tabulation of the exact solution. %
u = zeros(1,N);
for i = 1 : N
    u(i) = x(i)^4 + x(i)^2 - 1;
end

% Approximate solution. %
y = zeros(1,N);

%===================================%
% Matrix of the linear system.      %
%===================================%
system_matrix = zeros(N,N);

%===================================%
% Diagonals of the system matrix:   %
% adiag - lower diagonal            %
% bdiag - main diagonal             %
% cdiag - upper diagonal            %
%===================================%
adiag = zeros(1,N);
bdiag = zeros(1,N);
cdiag = zeros(1,N);

%===================================%
% Right side of the linear system.  %
%===================================%
right_side = zeros(N,1);

%===================================%
% Filling adiag with values.        %
%===================================%
adiag(1) = 0;
for i = 2 : N-1
    adiag(i) = k(x(i) - h/2) / (h^2);
end
adiag(end) = -k(x(N) - h/2) / h;

%===================================%
% Filling bdiag with values.        %
%===================================%
bdiag(1) = beta_1 + k(x(1) + h/2) / h + ...
    (h / 2) * q(x(1) + h/4);
for i = 2 : N-1
    bdiag(i) = -q(x(i)) - ...
                k(x(i) - h / 2) / (h^2) - ...
                k(x(i) + h / 2) / (h^2);
end
bdiag(end) = k(x(N) - h/2) / h + beta_2 + ...
    (h / 2) * q(x(N) - h/4);
%===================================%
% Filling cdiag with values.        %
%===================================%
cdiag(1) = -k(x(1) + h/2) / h;
for i = 2 : N-1
    cdiag(i) = k(x(i) + h/2) / (h^2);
end
cdiag(end) = 0;

%===================================%
% Filling right side with values.   %
%===================================%
right_side(1) = (h / 2) * f(x(1) + h/4) + mu_1;
for i = 2 : N-1
    right_side(i) = -f(x(i));
end
right_side(end) = (h / 2) * f(x(N) - h/4) + mu_2;

%==========================================%
% Filling the system matrix with values.   %
%==========================================%
for i0 = 1 : N
    for j0 = 1 : N
        if (i0 == j0 + 1)
            system_matrix(i0,j0) = adiag(i0);
        elseif (j0 == i0 + 1)
            system_matrix(i0,j0) = cdiag(i0);
        elseif(i0 == j0)
            system_matrix(i0,j0) = bdiag(i0);
        end
    end
end

%==============================%
% Thomas method (progonka)     %
%==============================%
y = Progon(system_matrix,right_side)';

%==============================%
% Absolute error.              %
%==============================%
abs_err = abs(u - y);
Max = max(max(abs_err));
disp(['Max. absolute error: ',num2str(Max)])

%==============================%
% Plot.                        %
%==============================%
figure(1)
% % %
subplot(1,2,1)
% % %
plot(x,y,'b','LineWidth',3)
hold on, grid on
plot(x,y,'g--','LineWidth',3)
set(gca,'FontSize',14)
xlabel('\bf{x}')
ylabel('\bf{y}')
legend('\it{Approximate solution}','\it{Exact solution}')
% % %
subplot(1,2,2)
% % %
plot(x,abs_err,'r:','LineWidth',3)
hold on, grid on
set(gca,'FontSize',14)
xlabel('\bf{x}')
ylabel('\bf{y}')
legend('\it{Absolute error}')

time = toc;
disp(['Elapsed time: ',num2str(time),' sec.']);
