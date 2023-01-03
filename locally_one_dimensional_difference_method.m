%===============================================%
% Nonlinear one-dimensional schema              %
% for the non-linear heat equation:             %
% u^2 * u_{t} = (u * u_{x})_{x} + u_{zz} + f    %
% 0 < x, z < 1, 0 <= t <= 1,                    %
% u(t=0) = sin(x) + cos(z),                     %
% u(x=0) = cos(z) + t^4,                        %
% u(x=1) = sin(1) + cos(z) + t^4,               %
% u(z=0) = sin(x) + 1 + t^4,                    %
% u(z=1) = sin(x) + cos(1) + t^4,               %
% f = 4 * t^3 * u^2 - cos(x)^2 + ...            %
%     u * sin(x) + cos(z)                       %
%===============================================%
% Exact solution for comparison:                %
% u(x,z,t) = sin(x) + cos(z) + t^4              %
%===============================================%
% Difference method: fully implicit locally     % 
% one-dimensional schema with truncation error  %
% O(hx^2 + hz^2 + tau).                         %
%===============================================%

clear; clc; close; tic;

%======================================%
% Introducing spatial and time grids.  %
%======================================%
ax = 0; bx = 1; hx = (bx - ax) / 60;
x = ax : hx : bx;
N = length(x);

az = 0; bz = 1; hz = (bz - az) / 50;
z = az : hz : bz;
M = length(z);

t0 = 0; tend = 1; tau = (tend - t0) / 600;
t = t0 : tau : tend;
T = length(t);

%===================================%
% Tabulating the exact solution.    %
%===================================%
ur = zeros(N,M,T);
for i = 1 : N
    for j = 1 : M
        for n = 1 : T
            ur(i,j,n) = sin(x(i)) + ...
                cos(z(j)) + t(n)^4;
        end
    end
end

%===================================
% Approximate solution 
%===================================
u = zeros(N,M,T);

%===================================
% Initial condition
%===================================
for i = 1 : N
    for j = 1 : M
        u(i,j,1) = sin(x(i)) + cos(z(j));
    end
end

%=====================%
% Locally 1D schema.  %
%=====================%
for n = 1 : T-1
    
    %===============================%
    % Solution on the layer n+1/2.  %
    %===============================%
    sol = zeros(N,M);
    
    %===============================%
    % Matrix of the linear system.  %
    %===============================%
    system_matrix = zeros(N,N);
    
    %================================%
    % 1D problem with respect to x.  %
    %================================%
    for j = 1 : M
        %========================%
        % Lower diagonal adiag.  %
        %========================%
        adiag = zeros(1,N);
        adiag(1) = 0;
        for ii = 2 : N-1
            adiag(ii) = -tau * 0.5 * (u(ii-1,j,n) + u(ii,j,n)) / ...
                        ((u(ii,j,n)^2) * (hx^2));
        end
        adiag(end) = 0;
        
        %========================%
        % Main diagonal bdiag.   %
        %========================%
        bdiag = zeros(1,N);
        bdiag(1) = 1;
        for ii = 2 : N-1
            bdiag(ii) = 1 + tau/((u(ii,j,n)^2) * (hx^2)) * ...
                        0.5 * (u(ii-1,j,n) + ...
                               2 * u(ii,j,n) + ...
                               u(ii+1,j,n));
        end
        bdiag(end) = 1;
        
        %========================%
        % Upper diagonal cdiag.  %
        %========================%
        cdiag = zeros(1,N);
        cdiag(1) = 0;
        cdiag(end) = 0;
        for ii = 2 : N-1
            cdiag(ii) = -tau * 0.5 * (u(ii+1,j,n) + u(ii,j,n)) / ...
                ((u(ii,j,n)^2) * (hx^2));
        end
        %===================%
        % Right-hand side.  %
        %===================%
        right_side = zeros(N,1);
        right_side(1) = cos(z(j)) + (t(n)+tau/2)^4;
        for ii = 2 : N-1
            right_side(ii) = u(ii,j,n) + tau/(2*(u(ii,j,n)^2)) * ...
                (4 * (t(n)+tau/2)^3 * u(ii,j,n)^2 - cos(x(ii))^2 + ...
                u(ii,j,n) * sin(x(ii)) + cos(z(j)));
        end
        right_side(end) = sin(1) + cos(z(j)) + (t(n)+tau/2)^4;
        
        %====================================%
        % Defining the matrix of the system. %
        %====================================%
        for i0 = 1 : N
            for j0 = 1 : N
                if(i0 == j0+1)
                    system_matrix(i0,j0) = adiag(i0);
                elseif(j0 == i0+1)
                    system_matrix(i0,j0) = cdiag(i0);
                elseif(i0 == j0)
                    system_matrix(i0,j0) = bdiag(i0);
                end
            end
        end
        
        %===========================================%
        % Progonka method for obtaining the values  %
        % of the solution on the layer n+1/2.       %
        %===========================================%
        sol(:,j) = Progon(system_matrix, right_side);
    end
    
    %===============================%
    % Matrix of the linear system.  %
    %===============================%
    system_matrix = zeros(M,M);
    
    %================================%
    % 1D problem with respect to z.  %
    %================================%
    for i = 1 : N
        %========================%
        % Lower diagonal adiag.  %
        %========================%
        adiag = zeros(1,M);
        adiag(1) = 0;
        for jj = 2 : M-1
            adiag(jj) = -tau / (((sol(i,jj))^2) * (hz^2));
        end
        adiag(end) = 0;
        
        %========================%
        % Main diagonal bdiag.   %
        %========================%
        bdiag = zeros(1,M);
        bdiag(1) = 1;
        for jj = 2 : M-1
            bdiag(jj) = 1 + 2 * tau/(((sol(i,jj))^2) * (hz^2));
        end
        bdiag(end) = 1;
        
        %========================%
        % Upper diagonal bdiag.  %
        %========================%
        cdiag = zeros(1,M);
        cdiag(1) = 0;
        for jj = 2 : M-1
            cdiag(jj) = -tau / (((sol(i,jj))^2) * (hz^2));
        end
        cdiag(end) = 0;
        
        %====================%
        % Right-hand side.   %
        %====================%
        right_side = zeros(M,1);
        right_side(1) = sin(x(i)) + 1 + (t(n+1))^4;
        for jj = 2 : M-1
            right_side(jj) = sol(i,jj) + tau/(2 * (sol(i,jj))^2) * ...
                (4 * (t(n)+tau/2)^3 * (sol(i,jj))^2 - cos(x(i))^2 + ...
                sol(i,jj) * sin(x(i)) + cos(z(jj)));
        end
        right_side(end) = sin(x(i)) + cos(1) + (t(n+1))^4;
        
        %====================================%
        % Defining the matrix of the system. %
        %====================================%
        for i0 = 1 : M
            for j0 = 1 : M
                if(i0 == j0+1)
                    system_matrix(i0,j0) = adiag(i0);
                elseif(j0 == i0+1)
                    system_matrix(i0,j0) = cdiag(i0);
                elseif(i0 == j0)
                    system_matrix(i0,j0) = bdiag(i0);
                end
            end
        end
        
        %=======================%
        % Progonka for u(n+1).  % 
        %=======================%
        u(i,:,n+1) = Progon(system_matrix, right_side);
    end
end

%===================%
% Absolute error.   %
%===================%
abs_err = abs(ur-u);

%========%
% Plot.  %
%========%
figure(1)
% % %
subplot(2,2,1)
% % %
plot(x,ur(:,end-3,end-10),'b','LineWidth',3)
hold on, grid on
plot(x,u(:,end-3,end-10),'g--','LineWidth',3)
set(gca,'FontSize',14)
xlabel(['\bf{x, z = }', num2str(z(end-3)), ...
    '\bf{, t = }', num2str(t(end-10))])
ylabel('\bf{u(x)}')
legend('\bf{Exact solution}','\bf{Approximate solution}')
title('\it{x-direction}')
% % %
subplot(2,2,2)
% % %
plot(x,abs_err(:,end-3,end-10),'r:','LineWidth',3)
hold on, grid on
set(gca,'FontSize',14)
xlabel(['\bf{x, z = }', num2str(z(end-3)), ...
    '\bf{, t = }', num2str(t(end-10))])
ylabel('\bf{Error}')
legend('\bf{|u(x) - y|}')
title('\it{x-direction}')
% % %
subplot(2,2,[3,4])
% % %
surf(z,x,u(:,:,end-10))
colorbar
set(gca,'FontSize',14)
xlabel('\it{x}')
ylabel('\it{z}')
zlabel('\it{u}')
legend('\bf{u(x,z,t)}')
title(['t = ', num2str(t(end-10))])

%=====================================%
% Printing important information.     %
%=====================================%
Max = max(max(max(abs_err)));
rel_err = max(max(max(abs_err)))/max(max(max(u)))*100;
display(['tau = ', num2str(tau)])
display(['Maximal absolute error: ',num2str(Max)])
display(['Relative error: ',num2str(rel_err)])
time = toc;
display(['Elapsed time: ', num2str(time)])
display('==================================')
