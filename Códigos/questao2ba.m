clear 
close all
clc

% Oscilador linear em espaços de fase
function ydot = osciladorL(t,y,zeta, Wn, gama,We)
ydot = [y(2);gama*sin(We*t)-2*zeta*Wn*y(2)-Wn^2*y(1)];
end

Wn = 2*pi;                      % frequência natural
zeta = 0.0;                     % fator de amortecimento 
gama = 2;                       % intensidade do forçamento
We = 1.2*Wn;                    % frequência de excitação
y0 = [0; 0];                    % condição inicial 

fun = @(t,y)osciladorL(t,y,zeta,Wn,gama,We);

tspan = [0, 40];                % intervalo de tempo da simulação
rtol = 1e-6;                    % tolerância relativa
atol = 1e-9;                    % tolerância absoluta
h0 = 1e-3;                      % passo inicial

stepfun = @rkdp45;

[T,Y,H] = rk45_adaptive(fun, tspan, y0, h0, rtol, atol, stepfun);  % função do passo adaptativo


% SOLUÇÃO NUMÉRICA
% Passo adaptativo conforme Hairer–Wanner:
function [T,Y,H] = rk45_adaptive(fun, tspan, y0, h0, rtol, atol, stepfun)

% Parâmetros do controlador
fac = 0.9;         % fator de segurança - o livro recomenda
facmin = 0.2;      % limite inferior do fator
facmax = 5.0;      % limite superior quando passo é ACEITO (esses limites evitam mudanças bruscas no passo)
q = 4;             % ordem menor do par de respostas (ordem do erro estimado)
maxSteps = 1e6;    % limite de steps para não ficar preso em um loop

% Parâmetros
t0 = tspan(1);                   % tempo inicial 
tf = tspan(2);                   % tempo final
t  = t0;
y  = y0(:);                      % vetor condições iniciais
n  = numel(y);                   % número de equações do sistema
 
% Vetorizar tolerâncias - viram um vetor com o valor imputado n vezes
if isscalar(rtol), rtol = rtol*ones(n,1); else, rtol = rtol(:); end
if isscalar(atol), atol = atol*ones(n,1); else, atol = atol(:); end

% armazenamento dinâmico
T = t;                     % vetor com os tempos aceitos
Y = y;                     % vetor com as soluções aceitas
H = [];                    % histórico dos passos utilizados

h = h0;                    % definição do passo 0

for k = 1:maxSteps
    if (t - tf) >= 0       % se o tempo passar o número de iterações máxima, para o loop
        break
    end
    
    if (t + h - tf) > 0    % se o próximo passo ultrapassar tf corrige h para coincidir com o tf
        h = tf - t;
    end

    % Um passo RK45
    [y4, y5] = stepfun(fun, h, t, y);

    % Erro escalado (Eq. 4.10 – 4.11)
    sc  = atol + max(abs(y), abs(y5)) .* rtol;     % fórmula livro (4.10)
    e   = y5 - y4;
    err = sqrt( mean( (e ./ sc).^2 ) );            % norma RMS - livro (4.11)

    if err <= 1 || err == 0      % ACEITA o passo
        t  = t + h;
        y  = y5;                 % avança no tempo e segue com a resposta de ordem 5
        T(end+1,1) = t;          % salva o tempo
        Y(:,end+1) = y;          % salva a solução
        H(end+1,1) = h;          % salva o passo

        % novo passo (Eq. 4.13) com facmax > 1 
        if err == 0
            facGrow = facmax;                        % fator crescido - cresce no máximo caso o erro seja zero
        else
            facGrow = fac * err^(-1/(q+1));          % ajusta o fator de crescimento máximo
        end

        h = h * min(facmax, max(facmin, facGrow));   % ajusta o novo passo

    else                                        % REJEITA o passo
        facRed = fac * err^(-1/(q+1));          % fator reduzido - repete o passo com h reduzido
        h = h * min(1.0, max(facmin, facRed));  % facmax = 1 no passo rejeitado
        % volta ao topo do laço sem avançar o tempo
    end
end

end

% Pontos de Poincaré
Te = 2*pi/We;                                                     % período
t_samp = T(1) + Te*(1:floor((T(end)-T(1))/Te));                   % criação de um vetor tempo que coincide com 'partes' do período de excitação

X_samp = interp1(T, Y(1,:), t_samp, 'pchip');                     % interpolação para pbter a resposta do deslocamento naquele tempo
V_samp = interp1(T, Y(2,:), t_samp, 'pchip');                     % interpolação para pbter a resposta da velocidade naquele tempo
poincare = [X_samp;
            V_samp];

% Regime permanente
rp = round(0.85*length(T));
rp_poincare = round(0.85*size(poincare,2));  % 85% do total de pontos de Poincaré

% Mapa de Poincaré
figure;
plot(Y(1,rp:end), Y(2,rp:end), 'b-', 'LineWidth', 1.2);                                      % Caminho contínuo
hold on;
plot(poincare(1,rp_poincare:end), poincare(2,rp_poincare:end), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5); % Pontos Poincaré
xlabel('$x(t)$', 'Interpreter', 'latex');
ylabel('$\dot{x}(t)$', 'Interpreter', 'latex');
title('Mapa de Poincaré por interpolação')

% n é tao preciso, existem outras formas de fazer que trazem menos erros,
% mas serviu para esse caso
