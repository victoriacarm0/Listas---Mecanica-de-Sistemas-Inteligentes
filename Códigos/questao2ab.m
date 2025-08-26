clear 
close all
clc

% Oscilador tipo Duffing em espaços de fase
function xdot = duffing(t,x,zeta,alfa,beta,gama,We)
xdot = [x(2);gama*sin(We*t)-zeta*x(2)+alfa*x(1)-beta*x(1)^3];
end

% Parâmetros do oscilador
alfa = -1.2;                       
beta = 0.3;
gama = 5;                      % amplitude de excitação
zeta = 0.05;                   % fator de amortecimento
We = 1;                      % frequência de excitação
Te = 2*pi/We;                  % período da excitação


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
y = 1;
for i = 1:length(tspan)-1 
    time = i*dt;
    xout = rk4(@(t,x)duffing(t,x,zeta,alfa,beta,gama,We), dt, time, xin);
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
xlabel('$x(t)$', 'Interpreter', 'latex');
ylabel('$\dot{x}(t)$', 'Interpreter', 'latex');
title('Mapa de Poincaré - Oscilador Duffing Forçado');
grid on;