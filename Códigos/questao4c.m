clc
clear
close all

% --- Parâmetros do sistema (da sua imagem 1) ---
z1 = 0.05;
z2 = 0.08;
a1 = -2;
a2 = -2;
b1 = 1;
b2 = 1.5;
rho = 0.5;
Os = 3;

% --- Função em espaço de estados ---
function dxdt = sistema2GDL(t,X,z1,z2,a1,a2,b1,b2,rho,Os)
    x1 = X(1); v1 = X(2);
    x2 = X(3); v2 = X(4);

    dx1 = v1;
    dv1 = -2*z1*v1 + 2*z2*(v2 - v1) - (1+a1)*x1 - b1*x1^3 + rho*Os^2*(x2 - x1);
    dx2 = v2;
    dv2 = (-(2*z2)*(v2 - v1) - a2*x2 - b2*x2^3 - rho*Os^2*(x2 - x1))/rho;

    dxdt = [dx1; dv1; dx2; dv2];
end

% Pontos de equilíbrio - encontrados na 3ª questão
eq_points = [ 1.0742  0   1.1109  0;
              0       0   0       0;
             -1.0742  0  -1.1109  0];

% Malha de condições iniciais 
x1_vals = linspace(-2,2,200);
x2_vals = linspace(-2,2,200);
[X1,X2] = meshgrid(x1_vals,x2_vals);

% Integração
tspan = [0 100];
tol = 1e-2; % tolerância de convergência
basin = zeros(size(X1));

for i = 1:numel(X1)
    % condição inicial: posição varia, velocidades = 0
    x0 = [X1(i); 0; X2(i); 0];
    [t,x] = ode45(@(t,X) sistema2GDL(t,X,z1,z2,a1,a2,b1,b2,rho,Os),tspan,x0);

    % estado final médio
    xf = mean(x(round(0.8*end):end,:),1);

    % distância aos equilíbrios
    dist = vecnorm(eq_points - xf,2,2);
    [dmin,idx] = min(dist);

    if dmin < tol
        basin(i) = idx;
    else
        basin(i) = 0; % sem convergência
    end
end

figure; hold on
basin_plot = reshape(basin,size(X1));

colors = [0.9 0.9 0.9;   % cinza = sem convergência
          0.3 0.6 1.0;   % azul = equilíbrio 1
          0.2 0.8 0.2;   % verde = equilíbrio 2
          1.0 0.2 0.2];  % vermelho = equilíbrio 3

image(x1_vals,x2_vals,basin_plot);
colormap(colors);
set(gca,'YDir','normal')
xlabel('$x_1(0)$','Interpreter','latex')
ylabel('$x_2(0)$','Interpreter','latex')
title('Bacia de Atração - Sistema 2 GDL')

% Marcar equilíbrios 
% extremos (estáveis) -> círculos
plot(eq_points([1 3],1), eq_points([1 3],3), 'ko', ...
    'MarkerFaceColor','y','MarkerSize',8)

% centro (instável) -> triângulo
plot(eq_points(2,1), eq_points(2,3), 'k^', ...
    'MarkerFaceColor','r','MarkerSize',9)

legend('Espirais Estáveis','Sela Instável')
hold off

% Bacia de Atração: x1 vs v1 
x1_vals = linspace(-2,2,200);
v1_vals = linspace(-5,5,200);
[X1,V1] = meshgrid(x1_vals, v1_vals);
basin_x1v1 = zeros(size(X1));

for i = 1:numel(X1)
    x0 = [X1(i); V1(i); 0; 0]; % x2 = 0, v2 = 0
    [t,x] = ode45(@(t,X) sistema2GDL(t,X,z1,z2,a1,a2,b1,b2,rho,Os),tspan,x0);
    xf = mean(x(round(0.8*end):end,:),1);
    dist = vecnorm(eq_points - xf,2,2);
    [dmin,idx] = min(dist);

    if dmin < tol
        basin_x1v1(i) = idx;
    else
        basin_x1v1(i) = 0;
    end
end

figure; hold on
image(x1_vals,v1_vals,reshape(basin_x1v1,size(X1)));
colormap(colors);
set(gca,'YDir','normal')
xlabel('$x_1(0)$','Interpreter','latex')
ylabel('$\dot{x}_1(0)$','Interpreter','latex')
title('Bacia de Atração - Sistema 2 GDL')

% Marcar equilíbrios 
% extremos (estáveis) -> círculos
plot(eq_points([1 3],1), eq_points([1 3],3), 'ko', ...
    'MarkerFaceColor','y','MarkerSize',8)

% centro (instável) -> triângulo
plot(eq_points(2,1), eq_points(2,3), 'k^', ...
    'MarkerFaceColor','r','MarkerSize',9)

legend('Espirais Estáveis','Sela Instável')

hold off

% Bacia de Atração: x2 vs v2 
x2_vals = linspace(-2,2,200);
v2_vals = linspace(-5,5,200);
[X2,V2] = meshgrid(x2_vals, v2_vals);
basin_x2v2 = zeros(size(X2));

for i = 1:numel(X2)
    x0 = [0; 0; X2(i); V2(i)]; % x1 = 0, v1 = 0
    [t,x] = ode45(@(t,X) sistema2GDL(t,X,z1,z2,a1,a2,b1,b2,rho,Os),tspan,x0);
    xf = mean(x(round(0.8*end):end,:),1);
    dist = vecnorm(eq_points - xf,2,2);
    [dmin,idx] = min(dist);

    if dmin < tol
        basin_x2v2(i) = idx;
    else
        basin_x2v2(i) = 0;
    end
end

figure; hold on
image(x2_vals,v2_vals,reshape(basin_x2v2,size(X2)));
colormap(colors);
set(gca,'YDir','normal')
xlabel('$x_2$','Interpreter','latex')
ylabel('$\dot{x}_2$','Interpreter','latex')
title('Bacia de Atração - Sistema 2 GDL')

% --- Marcar equilíbrios ---
% extremos (estáveis) -> círculos
plot(eq_points([1 3],1), eq_points([1 3],3), 'ko', ...
    'MarkerFaceColor','y','MarkerSize',8)

% centro (instável) -> triângulo
plot(eq_points(2,1), eq_points(2,3), 'k^', ...
    'MarkerFaceColor','r','MarkerSize',9)

legend('Espirais Estáveis','Sela Instável')

hold off