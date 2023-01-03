%===================================%
% Using ode45() built-in funtion    %
% to solve an ODE of 3-rd order by  %
% reducing to a 1-st order system   %
% of ODEs.                          %
%===================================%
% Initial value problem:            %
% y''' + y = sin(x), 0 < x < 5,     %
% y(0) = 0, y'(0) = 1, y''(0) = 1   %
%===================================%
% The exact solution is obtained    %
% using dsolve() built-in function. %
%===================================%

clear; clc;

%====================================%
% Right-hand side of the ODE system. % 
%====================================%
fun = @(x,u) [u(2); u(3); sin(x)-u(1)];

[xGrid, yGrid] = ode45(fun, [0, 5], [0; 1; 1]);
hold on 
grid on
axis([0, 5, 0, 6])
set(gca, 'FontName', 'Times', 'FontSize', 14)
title('$$ y'''''' + y = \sin x, \quad y(0)=0, \quad y''(0)=1, \quad y''''(0)=1 $$', ...
    'interpreter', 'latex')
xlabel('$$ x \in [0; \, 5] $$', 'interpreter', 'latex')
ylabel('$$ y \in [0; \, 6] $$', 'interpreter', 'latex')

%=============================%
% Exact solving by dsolve().  %
%=============================%
uExact = dsolve('D3y + y = sin(x), y(0)=0, Dy(0)=1, D2y(0)=1', 'x');
uExact = sym(uExact); % symbolic in order to use simplify() and pretty()
uExact = simplify(uExact);

%==============================%
% Printing the exact solution  %
% in explicit form.            %
%==============================%
disp('Exact solution: y(x)=')
pretty(uExact)

%===========%
% Plot.     %
%===========%
uGrid = subs(uExact, 'x', xGrid);
plot(xGrid, yGrid(:, 1), 'b', 'LineWidth', 3)
plot(xGrid,uGrid,'c:','LineWidth',3)
legend('\it{Numerical solution}', '\it{Exact solution}')
