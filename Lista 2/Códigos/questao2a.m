% brinson_SMA_cycle_corrigido.m
% Ciclo memória de forma 
clear; close all; clc;

% Parâmetros
Ca = 67e3;    % MPa
Cm = 26.3e3;  % MPa
Theta = 0.55; % MPa/°C
eps_R = 0.067;

Mf = 9.0;  Ms = 18.4;
As = 34.5; Af = 49.0;

C_M = 8.0;   % MPa/°C
C_A = 13.8;  % MPa/°C
sigma_s_cr = 100; % MPa
sigma_f_cr = 170; % MPa

%% Tempo e carregamento
Tend = 8; dt = 0.02;
time = (0:dt:Tend)';

% sigma: aplica 0→200 MPa (0–2s), tira (2–4s), depois zero
sigma_max = 200; % MPa
sigma = zeros(size(time));
for k = 1:length(time)
    t = time(k);
    if t<=2
        sigma(k) = sigma_max*(t/2);
    elseif t<=4
        sigma(k) = sigma_max*(1-(t-2)/2);
    else
        sigma(k) = 0;
    end
end

% temperatura: aplica 0→50 °C (4–6s), tira (6–8s), senão zero
theta = zeros(size(time));
for k = 1:length(time)
    t = time(k);
    if t<=4
        theta(k) = 0;
    elseif t<=6
        theta(k) = 50*(t-4)/2;
    else
        theta(k) = 50*(1-(t-6)/2);
    end
end

%% Estados
eps = zeros(size(time));
betaS = zeros(size(time));
betaT = zeros(size(time));

% Condição inicial: martensita maclada
betaS(1) = 0; 
betaT(1) = 1;
eps(1) = 0;

sigma0 = 0; 
eps0 = 0; 
theta0 = 0;
betaS0 = betaS(1); 
betaT0 = betaT(1);
beta0  = betaS0 + betaT0;

%% Loop de integração
for k = 2:length(time)
    s_app = sigma(k);
    th = theta(k);

    % Atualiza frações de fase com base no estado anterior
    [bS_new, bT_new] = update_betas(betaS(k-1), betaT(k-1), s_app, th, Mf, Ms, As, Af, sigma_s_cr, sigma_f_cr, C_M, C_A);
    betaS(k) = bS_new;
    betaT(k) = bT_new;
    beta = betaS(k) + betaT(k);

    % Constante elástica instantânea
    Cbeta = Ca + beta*(Cm - Ca);
    Cbeta0 = Ca + beta0*(Cm - Ca);

    % Gamas
    gammaS = - Cbeta*eps_R;
    gammaT = - Cbeta0*eps_R;

    % Equação constitutiva corrigida (inclui deformação elástica sempre)
    eps(k) = (Cbeta0*eps0 - gammaS*betaS(k) + gammaT*betaS0 - Theta*(th - theta0) + (s_app - sigma0))/Cbeta;
       
    % Atualização dos betas "iniciais da transformação"
    if th >= Af && s_app < sigma_s_cr                     % austenita pura
        betaS0 = 0; betaT0 = 0;
        betaS(k) = 0; betaT(k) = 0;
    elseif th >= Af && s_app >= sigma_f_cr                 % martensita desmacleada
        betaS0 = 1; betaT0 = 0;
        betaS(k) = 1; betaT(k) = 0;
    elseif th <= As && s_app <= sigma_s_cr && sigma(k-1) <= sigma_s_cr  % martensita maclada
        betaS0 = 0; betaT0 = 1;
        %betaS(k) = 0; betaT(k) = 1;
    %elseif th <= As && s_app <= sigma_s_cr %&& sigma(k-1) >= sigma_s_cr
        %betaS0 = 1; betaT0 = 0;
    end
end

% Gráficos
figure
%subplot(2,2,1);
yyaxis left; plot(time,sigma,'k','LineWidth',1.2); ylabel('\sigma (MPa)');
yyaxis right; plot(time,theta,'r','LineWidth',1.2); ylabel('T (°C)');
xlabel('t [s]'); grid on; title('Carregamento aplicado');

figure;
%subplot(2,2,2);
plot(time,betaS,'m','LineWidth',1.5); hold on;
plot(time,betaT,'y','LineWidth',1.5);
xlabel('t [s]'); ylabel('\beta'); legend('\beta_S','\beta_T');
title('Frações de fase'); grid on;

figure;
%subplot(2,2,3);
plot(eps,sigma,'b','LineWidth',1.5);
xlabel('\epsilon'); ylabel('\sigma (MPa)');
title('Laço pseudoelástico/memória de forma'); grid on;

figure;
%subplot(2,2,4);
plot3(theta,eps,sigma,'b','LineWidth',1.2);
xlabel('T [°C]'); ylabel('\epsilon'); zlabel('\sigma (MPa)');
ylim([-0.02 0.08]);
title('Trajetória 3D'); grid on; view(45,20);

%% Função auxiliar corrigida
function [betaS_new, betaT_new] = update_betas(betaS_prev, betaT_prev, sigma_eff, theta, Mf, Ms, As, Af, sigma_s_cr, sigma_f_cr, C_M, C_A)

    betaS_new = betaS_prev; 
    betaT_new = betaT_prev; 
    beta_prev = betaS_prev + betaT_prev;

    % Martensita desmacleada
    if (theta > Ms) && (sigma_eff > (sigma_s_cr + C_M*(theta - Ms))) && (sigma_eff < (sigma_f_cr + C_M*(theta - Ms)))
        arg = pi/(sigma_s_cr - sigma_f_cr) * (sigma_eff - sigma_f_cr - C_M*(theta - Ms));
        betaS_new = (1-betaS_prev)/2*cos(arg) + (1+betaS_prev)/2;
        betaT_new = betaT_prev - betaT_prev/(1-betaS_prev)*(betaS_new-betaS_prev);
    end
    
    if (theta < Ms) && ( sigma_eff > sigma_s_cr ) && ( sigma_eff < sigma_f_cr )
        denom = (sigma_s_cr - sigma_f_cr);
        arg = pi/denom * ( sigma_eff - sigma_f_cr );
        betaS_new = (1 - betaS_prev)/2 * cos(arg) + (1 + betaS_prev)/2;
        betaT_new = betaT_prev - betaT_prev/(1 - betaS_prev) * (betaS_new - betaS_prev);

        % Delta_Tbeta se Mf < theta < Ms
        if (theta > Mf) && (theta < Ms)
            DTbeta = (1 - betaS_prev)/2 * ( cos( pi/(Ms - Mf) * (theta - Mf) ) + 1 );
            betaT_new = betaT_new + DTbeta;
        end
    end

    % Austenita
    if (theta > As) && (sigma_eff > (C_A*(theta - Af))) && (sigma_eff < (C_A*(theta - As))) && (theta < Af)
        arg = pi/(Af-As)*(theta - As - sigma_eff/C_A);
        beta_new = beta_prev/2*(cos(arg)+1);
        betaS_new = betaS_prev - betaS_prev/beta_prev*(beta_prev-beta_new);
        betaT_new = beta_new - betaS_new;
    end

    % Limitadores
    betaS_new = min(max(betaS_new,0),1);
    betaT_new = min(max(betaT_new,0),1);
end