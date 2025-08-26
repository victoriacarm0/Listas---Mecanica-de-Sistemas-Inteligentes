clc
clear
close all

% Função do sistema do Pêndulo Amortecido
function dxdt = pendulo(t,x,zeta)
    dxdt = [x(2); -2*zeta*x(2) - sin(x(1))];
end

% Parâmetro de amortecimento
zeta = 0.1;

% Pontos de equilíbrio
eq_points = [0 0; pi 0; -pi 0];   % (0,0) estável, (±pi,0) instável

% Malha de condições iniciais
phi_vals  = linspace(-3, 3, 100);   % ângulo inicial
dphi_vals = linspace(-3, 3, 100);   % velocidade inicial
[X1, X2] = meshgrid(phi_vals, dphi_vals);

% Parâmetros de integração
tspan = [0 200];   % tempo de simulação
tol = 1;       % tolerância para identificar convergência

% Matriz para armazenar bacia
basin = zeros(size(X1));

% Loop sobre condições iniciais x1 x x2
for i = 1:numel(X1)
    x0 = [X1(i); X2(i)];
    [t, x] = ode45(@(t,x) pendulo(t,x,zeta), tspan, x0);
    
    xf = x(end,:);   % estado final
    
    % Distância até pontos de equilíbrio
    dist = vecnorm(eq_points - xf, 2, 2);
    [dmin, idx] = min(dist);
    
    if dmin < tol
        basin(i) = idx;   % atraiu para esse equilíbrio
    else
        basin(i) = 0;     % não atraiu para nenhum
    end
end


% Plotagem
figure; hold on
basin_plot = reshape(basin, size(X1));

% Cores: 0 = cinza, (0,0) azul, (pi,0) vermelho, (-pi,0) vermelho
colors = [0.9 0.9 0.9;   % sem convergência
          0.6 0.8 1.0;   % azul claro = estável
          1.0 0.6 0.6;   % vermelho claro = instável
          1.0 0.6 0.6];  % idem para -pi

image(phi_vals, dphi_vals, basin_plot);
colormap(colors);
set(gca,'YDir','normal')
xlabel('$\varphi$ (rad)', 'Interpreter', 'latex');
ylabel('$\dot{\varphi}$ (rad/s)', 'Interpreter', 'latex');
title('Bacia de Atração do Pêndulo Amortecido')

% Marcar pontos de equilíbrio
plot(0,0,'ro','MarkerFaceColor','b','MarkerSize',8)   % Estável
plot(pi,0,'rs','MarkerFaceColor','r','MarkerSize',8)  % Instável
plot(-pi,0,'rs','MarkerFaceColor','r','MarkerSize',8) % Instável

legend('Estável (0,0)','Instável (±π,0)')

hold off
