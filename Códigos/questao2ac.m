clear 
close all
clc

% Pêndulo simples em espaços de fase
function xdot = pendulo(t,x, Wn, gama, zeta, We)
xdot = [x(2);gama*sin(We*t)-Wn^2*sin(x(1))-zeta*x(2)];
end

% Parâmetros do oscilador
gama = 1;                       % amplitude de excitação
zeta = 0.05;                    % fator de amortecimento
Wn = 2*pi;                      % frequência natural
We = 2*Wn;                      % frequência de excitação
Te = 2*pi/We;                   % período da excitação


% Parâmetros para simulação 
x0 = [0; 0];                   % condição inicial x e xponto
dt = 0.01;                     % passo desejado 
n = round(Te/dt);
dt = Te / n;                   % passo ajustado
tf = 100*Te;
tspan = 0 : dt : tf;           % tempo de simulação

X = [];
X(:,1) = x0;
xin = x0;
poincare = [];
y = 2;
for i = 1:length(tspan)-1 
    time = i*dt;
    xout = rk4(@(t,x)pendulo(t,x, Wn, gama, zeta, We), dt, time, xin);
    X = [X xout];
    
    if i == y 
        poincare(:, end+1) = xout;
        y = i + n;
    end

    xin = xout;
end

% Regime permanente
rp = round(0.85*length(tspan));
rp_poincare = round(0.85*size(poincare,2));  % 85% do total de pontos de Poincaré

% Plotando
figure;
plot(X(1,rp:end), X(2,rp:end), 'b-', 'LineWidth', 1.2);                                       % Caminho contínuo
hold on;
plot(poincare(1,rp_poincare:end), poincare(2, rp_poincare:end), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5); % Pontos Poincaré
xlabel('$\theta$', 'Interpreter', 'latex');
ylabel('$\dot{\theta}$', 'Interpreter', 'latex');
title('Mapa de Poincaré - Pêndulo Forçado');
grid on;
figure;
plot(tspan, X(1,:), 'm', 'LineWidth', 0.75);