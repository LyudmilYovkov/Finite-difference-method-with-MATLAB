%=========================================================%
% u_{t} = (u^sigma * u_{x})_{x},                          %
% sigma > 0, 0 < x < X, t > 0,                            %
% IC: u(x,0) = 0, 0 <= x <= X                             %
% left BC: u(0,t) = t^(1/sigma), t >= 0                   %
% right BC: u(X,t) = 0, t >= 0                            %
%=========================================================%
% Exact solution: u(x,t)                                  %
% if (0 <= x <= D * t)                                    %
%   u(x,t) = (1 / D)^(1/sigma) * (D * t - x)^(1/sigma)    %
% else if (x > D * t)                                     %
%   u(x,t) = 0                                            %
% end                                                     %
%=========================================================%
% Constants: sigma, D, X                                  %
% sigma = [0.5; 1; 2]                                     %
% D = 1 / sqrt(sigma)                                     %  
% X - positive real number                                %
%=========================================================%
% Difference method: fully implicit scheme                %
%=========================================================%

clear; clc;

% Parameters.
sigma = 0.5;
D = 1 / sqrt(sigma);
X = 4;

% Heat coefficient.
k = @(u) (u).^(sigma);

% x-grid.
x0 = 0; 
x = linspace(x0,X,300);
h = x(2) - x(1);
N = numel(x);

% t-grid.
t0 = 0;
T = 1;
t = linspace(t0,T,640);
tau = t(2) - t(1);
M = numel(t);

% Tabulation of the exact solution.
u = zeros(N,M);
for i = 1 : N
    for j = 1 : M
        if (0 <= x(i) && x(i) <= D * t(j))
            u(i,j) = (1 / D)^(1/sigma) * ...
                     (D * t(j) - x(i))^(1/sigma);
        elseif(x(i) > D * t(j))
            u(i,j) = 0;
        end
    end
end

% Approximate solution.
y = zeros(N,M);

% Initial condition.
for i = 1 : N
    y(i,1) = 0; 
end

% Boundary conditions.
for j = 1 : M
    y(1,j) = (t(j))^(1/sigma); 
    y(N,j) = 0;
end

% Matrix of the system.
system_matrix = zeros(N,N);

% Right side.
right_side = zeros(N,1);

%==============================%
% Diagonal vectors:            %
% adiag - lower diagonal       %
% bdiag - main diagonal        %
% cdiag - upper diagonal       %
%==============================%
adiag = zeros(1,N);
bdiag = zeros(1,N);
cdiag = zeros(1,N);

% Finite difference method.
for j = 1 : M-1     
    %============================%
    % Filling adiag with values. %
    %============================%
    adiag(1) = 0;
    for i = 2 : N-1
        adiag(i) = ...
            -k( ( y(i-1,j) + y(i,j) ) / 2 ) / (h^2);
    end
    adiag(end) = 0;
    
    %============================%
    % Filling bdiag with values. %
    %============================%
    bdiag(1) = 1;
    for i = 2 : N-1
        bdiag(i) = 1/tau + ...
            k( ( y(i-1,j) + y(i,j) ) / 2 ) / (h^2) + ...
            k( ( y(i,j) + y(i+1,j) ) / 2 ) / (h^2);
    end
    bdiag(end) = 1;
    
    %============================%
    % Filling cdiag with values. %
    %============================%
    cdiag(1) = 0;
    for i = 2 : N-1
        cdiag(i) = ...
            -k( ( y(i,j) + y(i+1,j) ) / 2 ) / (h^2);
    end
    cdiag(end) = 0;
    
    %============================%
    % Filling the right side of  %
    % the system with values.    %
    %============================%
    right_side = zeros(N,1);
    right_side(1) = (t(j))^(1/sigma);
    for i = 2 : N-1
        right_side(i) = (1/tau) * y(i,j);
    end
    right_side(end) = 0;
    
    %============================%
    % Filling the matrix of the  %
    % system with values.        %
    %============================%
    system_matrix = zeros(N,N);
    for i0 = 1 : N
        for j0 = 1 : N
            if(i0 == j0 + 1)
                system_matrix(i0,j0) = adiag(i0);
            elseif(j0 == i0 + 1)
                system_matrix(i0,j0) = cdiag(i0);
            elseif(i0 == j0)
                system_matrix(i0,j0) = bdiag(i0);
            end
        end
    end
        
    %============================%
    % Using Thomas method for    %
    % finding y(j+1).            %
    %============================%
    y(:,j+1) = ...
        Progon(system_matrix,right_side);
end

% Absolute Error.
abs_err = abs(u - y);
Max = max(max(abs_err));
disp(['Max. absolute error: ',num2str(Max)])

% Plot at a given time.
figure(1)
time_layer = round(M/2);
% % %
subplot(1,2,1)
% % %
plot(x,y(:,time_layer),'b','LineWidth',3)
hold on, grid on
plot(x,u(:,time_layer),'g--','LineWidth',3)
set(gca,'FontSize',14)
xlabel('\bf{x}')
ylabel('\bf{y}')
legend('\it{Approximate solution}','\it{Exact solution}')
title('$  \sigma = 0.5 $','interpreter','latex')
% % %
subplot(1,2,2)
% % %
plot(x,abs_err(:,time_layer),'r:','LineWidth',3)
hold on, grid on
set(gca,'FontSize',14)
xlabel('\bf{x}')
ylabel('\bf{y}')
legend('\it{Absolute error}')
title('$  \sigma = 0.5 $','interpreter','latex')