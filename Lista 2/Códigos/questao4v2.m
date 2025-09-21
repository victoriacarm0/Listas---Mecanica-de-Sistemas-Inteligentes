clear; clc; close all;

% Parâmetros fixos
zeta = 0.025;        % Amortecimento adimensional
A = 5;               % Força de excitação (2.5 -> 9.81)
a = 15;              % [N/(K.m.kg)]
b = 60e4;            % [N/(m^3.kg)]
Ta = 313;            % [K]
Tm = 287;            % [K]
T  = 300;            % Temperatura atual [K]

% Frequências de excitação para FRF
omega_vec = linspace(1,100,50);   % 50 pontos de varredura

% Vetores para guardar amplitudes
amp_SMA = zeros(size(omega_vec));
amp_lin = zeros(size(omega_vec));

% Loop sobre frequências
for i = 1:length(omega_vec)
    omega = omega_vec(i);

    % Tempo de simulação (mais longo para deixar atingir regime)
    tspan = linspace(0,50,5000);
    x0 = [0 0];

    % Integra SMA
    [t, y_SMA] = ode45(@(t,y) SMA_oscillator(t,y,zeta,a,b,T,Ta,Tm,A,omega), tspan, x0);
    x_SMA = y_SMA(:,1);

    % Integra Linear
    [~, y_lin] = ode45(@(t,y) linear_oscillator(t,y,zeta,a,T,Tm,A,omega), tspan, x0);
    x_lin = y_lin(:,1);

    % Pega só a parte final (regime permanente)
    idx = round(length(t)*0.8):length(t);
    amp_SMA(i) = max(abs(x_SMA(idx)));
    amp_lin(i) = max(abs(x_lin(idx)));
end

% Plot FRF
figure;
plot(omega_vec, amp_SMA,'b','LineWidth',1.5); hold on;
plot(omega_vec, amp_lin,'r--','LineWidth',1.5);
xlabel('\omega [rad/s]');
ylabel('Amplitude máxima |x(t)| [m]');
title('Varredura em frequência: Oscilador SMA vs Linear');
legend('Oscilador SMA','Oscilador Linear');
grid on;

% Plot resposta no tempo - freq específica
omega_test = 14;             % escolhe a frequência [rad/s]
tspan = linspace(0,20,2000);
x0 = [0 0];

% Simulação SMA
[t, y_SMA] = ode45(@(t,y) SMA_oscillator(t,y,zeta,a,b,T,Ta,Tm,A,omega_test), tspan, x0);
x_SMA = y_SMA(:,1);

% Simulação Linear
[~, y_lin] = ode45(@(t,y) linear_oscillator(t,y,zeta,a,T,Tm,A,omega_test), tspan, x0);
x_lin = y_lin(:,1);

% Plot das respostas
figure;
plot(t,x_SMA,'b','LineWidth',1.5); hold on;
plot(t,x_lin,'r--','LineWidth',1.5);
xlabel('Tempo [s]');
ylabel('x(t) [m]');
title(['Resposta no tempo para \omega = ' num2str(omega_test) ' rad/s']);
legend('Oscilador SMA','Oscilador Linear');
grid on;

% Funções dos sistemas
function dydt = SMA_oscillator(t,y,zeta,a,b,T,Ta,Tm,A,omega)
    x = y(1); dx = y(2);
    k1 = a*(T - Tm);                
    k3 = -b;
    k5 = (b^2)/(4*a*(Ta - Tm));
    ddx = -2*zeta*dx - k1*x + k3*x^3 - k5*x^5 + A*sin(omega*t);
    dydt = [dx; ddx];
end

function dydt = linear_oscillator(t,y,zeta,a,T,Tm,A,omega)
    x = y(1);
    dx = y(2);
    k_lin = a*(T - Tm);
    ddx = -2*zeta*sqrt(k_lin)*dx - k_lin*x + A*sin(omega*t);
    dydt = [dx; ddx];
end
