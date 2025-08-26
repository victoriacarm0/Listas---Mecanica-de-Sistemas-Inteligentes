clear 
close all
clc

% Oscilador linear em espaços de fase
function xdot = osciladorL(t,x,zeta, Wn, gama,We)
xdot = [x(2);gama*sin(We*t)-2*zeta*Wn*x(2)-Wn^2*x(1)];
end

% Parâmetros do oscilador
zeta = 2;                    % amortecimento
gama = 2;                       % amplitude de excitação
Wn = 2*pi;                      % frequência natural
We = 1.2*Wn;                     % frequência de excitação
Wd = Wn*sqrt(1-zeta^2);         % frequência natural amortecida

% Parâmetros para simulação 
dt_list = [0.1, 0.05, 0.01]; 
cores = ['b', 'r', 'g'];       % cores para diferenciar

figure;
subplot(2,1,1); hold on;    % para os gráficos virem em uma página só
subplot(2,1,2); hold on;

% Variação de dt e aplicação do rk4
for j = 1:length(dt_list)
    dt = dt_list(j);
    tspan = 0:dt:30;
 
    % condições iniciais
    x0 = [1; 2]; 
    X = x0;
    xin = x0;

    % integração numérica (RK4)
    for i = 1:(length(tspan)-1)
        time = tspan(i);
        xout = rk4(@(t,x)osciladorL(t,x,zeta,Wn,gama,We), dt, time, xin);
        X = [X xout];
        xin = xout;
    end

   % solução analítica - solução particular + homogênea
   % solução particular 
    xp_amp = gama / sqrt( (Wn^2 - We^2)^2 + (2*zeta*Wn*We)^2 );
    phi = atan2( 2*zeta*Wn*We, (Wn^2 - We^2) );
    xp = xp_amp * sin(We*tspan - phi);

   % solução homogênea 
    if zeta == 0
    % Caso não amortecido
     C1 = x0(1) - xp(1);
     C2 = (x0(2) - We*xp_amp*cos(phi)) / Wn;
     xh = C1*cos(Wn*tspan) + C2*sin(Wn*tspan);

    elseif zeta < 1
    % Subamortecido
     Wd = Wn*sqrt(1 - zeta^2);
     C1 = x0(1) - xp(1);
     C2 = (x0(2) - We*xp_amp*cos(phi) + zeta*Wn*C1) / Wd;
     xh = exp(-zeta*Wn*tspan) .* ( C1*cos(Wd*tspan) + C2*sin(Wd*tspan) );

    elseif abs(zeta - 1) < 1e-6
    % Criticamente amortecido
     C1 = x0(1) - xp(1);
     C2 = x0(2) - We*xp_amp*cos(phi) + Wn*C1;
     xh = (C1 + C2*tspan) .* exp(-Wn*tspan);

    elseif zeta > 1
    % Superamortecido
     s1 = -Wn*(zeta - sqrt(zeta^2 - 1));
     s2 = -Wn*(zeta + sqrt(zeta^2 - 1));
     A = [1, 1; s1, s2] \ ([x0(1) - xp(1); x0(2) - We*xp_amp*cos(phi)]);
     C1 = A(1); C2 = A(2);
     xh = C1*exp(s1*tspan) + C2*exp(s2*tspan);
    end

    x_analit = xh + xp;

    % Erro absoluto
    erro = abs(X(1,:) - x_analit);

    % Plota cada dt com cor diferente
    subplot(2,1,1);
    plot(tspan, X(1,:), cores(j), 'DisplayName', ['Num dt=' num2str(dt)]);

    subplot(2,1,2);
    plot(tspan, erro, cores(j), 'DisplayName', ['Erro dt=' num2str(dt)]);
end

% Formatação dos gráficos
subplot(2,1,1); grid on; legend; 
xlabel('Tempo (s)');
ylabel('Deslocamento x(t)');
title('Comparação entre solução numérica e analítica');
plot(tspan, x_analit, '--m', 'LineWidth', 0.75);

subplot(2,1,2); grid on; legend; title('Erro absoluto');
xlabel('Tempo (s)');
ylabel('Erro absoluto');
title('Erro Absoluto');
grid on;