clc
clear
close all

% --- Função do sistema de Duffing ---
function xout = duffing(t,x,alfa,beta,zeta,gama,We) 
    xout = [x(2) ; -2*zeta*x(2) + alfa*x(1) - beta*x(1)^3 + gama*sin(We*t)];
end

% Parâmetros do sistema
alfa = 1;
beta = 1;
zeta = 0.1;            % amortecimento
gama = 0;              % sem excitação externa
We = 1;                % não será usado pois gamma=0

% Pontos de equilíbrio (x1, x2)
eq_points = [0               0; 
            sqrt(alfa/beta)  0; 
            -sqrt(alfa/beta) 0];

% Malha de condições iniciais
x1_vals = linspace(-2, 2, 200);         % variação de posição inicial
x2_vals = linspace(-2, 2, 200);         % variação de velocidade inicial
[X1, X2] = meshgrid(x1_vals, x2_vals);  % cria duas matrizes 2D que representam todos os pares possíveis x e v

% Parâmetros de integração
tspan = [0 100];  % tempo grande para garantir convergência
tol = 1e-3;       % tolerância para identificar equilíbrio

% Matriz para armazenar o índice do ponto de equilíbrio
basin = zeros(size(X1));

% Loop sobre condições iniciais
for i = 1:numel(X1)
    x0 = [X1(i); X2(i)];
    [t, x] = ode45(@(t,x) duffing(t,x,alfa,beta,zeta,gama,We), tspan, x0);
    
    % Pega estado final
    xf = x(end,:);
    
    % Verifica para qual equilíbrio foi
    dist = vecnorm(eq_points - xf, 2, 2);   % distância Euclidiana entre o estado final e cada equilíbrio 
    [~, idx] = min(dist);                   % pega a menor distância e o índice
    
    % Se distância pequena, considera convergência
    if dist(idx) < tol                      % verifica se realmente tá bem perto, se sim, convergiu, se não, não
        basin(i) = idx;
    else
        basin(i) = 0;                       % não convergiu a nenhum 
    end
end

% Plotagem
figure; hold on
basin_plot = reshape(basin, size(X1));

% Definindo cores mais claras para as bacias
colors = [0.9 0.9 0.9;   % cinza claro = não convergiu
          1.0 0.6 0.6;   % vermelho claro = (0,0)
          0.6 0.8 1.0;   % azul claro = (sqrt(alpha/beta),0)
          0.6 1.0 0.6];  % verde claro = (-sqrt(alpha/beta),0)

image(x1_vals, x2_vals, basin_plot);
colormap(colors);
set(gca,'YDir','normal')
xlabel('$x(t)$', 'Interpreter', 'latex');
ylabel('$\dot{x}(t)$', 'Interpreter', 'latex');
title('Bacia de Atração - Oscilador de Duffing')

% Marcar pontos de equilíbrio
plot(eq_points(1,1), eq_points(1,2), 'k^','MarkerFaceColor','y','MarkerSize',8) % Sela (instável)
plot(eq_points(2,1), eq_points(2,2), 'ko','MarkerFaceColor','b','MarkerSize',8) % Estável
plot(eq_points(3,1), eq_points(3,2), 'ko','MarkerFaceColor','b','MarkerSize',8) % Estável

legend('Sela (instável)','Espirais estáveis')
hold off

