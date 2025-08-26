clear 
close all
clc

% Oscilador linear em espaços de fase
function ydot = osciladorL(t,y,zeta, Wn, gama,We)
ydot = [y(2);gama*sin(We*t)-2*zeta*Wn*y(2)-Wn^2*y(1)];
end

Wn = 2*pi;                      % frequência natural
zeta = 0.05;                    % fator de amortecimento 
gama = 1;                       % intensidade do forçamento
We = 1.2*Wn;                    % frequência de excitação
Wd = Wn*sqrt(1-zeta^2);         % frequência natural amortecida
y0 = [1; 2];                    % condição inicial 

fun = @(t,y)osciladorL(t,y,zeta,Wn,gama,We);

tspan = [0, 100];                % intervalo de tempo da simulação
rtol = 1e-6;                    % tolerância relativa
atol = 1e-9;                    % tolerância absoluta
h0 = 1e-3;                      % passo inicial

stepfun = @rkdp45;

[T,Y,H,stats] = rk45_adaptive(fun, tspan, y0, h0, rtol, atol, stepfun);  % função do passo adaptativo

fprintf('Passos aceitos: %d | rejeitados: %d\n', stats.accepted, stats.rejected);  % plota n° de passos aceitos/rejeitados

figure;
subplot(2,1,1); hold on;    % para os gráficos virem em uma página só
subplot(2,1,2); hold on;

% SOLUÇÃO NUMÉRICA
% Passo adaptativo conforme Hairer–Wanner:
function [T,Y,H,stats] = rk45_adaptive(fun, tspan, y0, h0, rtol, atol, stepfun)

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
nacc = 0;                  % contador de passos aceitos
nrej = 0;                  % contador de passos rejeitados

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
        nacc = nacc + 1;         % adiciona 1 ao contador de passos aceitos

        % novo passo (Eq. 4.13) com facmax > 1 
        if err == 0
            facGrow = facmax;                        % fator crescido - cresce no máximo caso o erro seja zero
        else
            facGrow = fac * err^(-1/(q+1));          % ajusta o fator de crescimento máximo
        end

        h = h * min(facmax, max(facmin, facGrow));   % ajusta o novo passo

    else                                        % REJEITA o passo
        nrej = nrej + 1;                        % adiciona 1 ao contador de passos rejeitados
        facRed = fac * err^(-1/(q+1));          % fator reduzido - repete o passo com h reduzido
        h = h * min(1.0, max(facmin, facRed));  % facmax = 1 no passo rejeitado
        % volta ao topo do laço sem avançar o tempo
    end
end

stats.accepted  = nacc;
stats.rejected  = nrej;
end

% SOLUÇÃO ANALÍTICA
% solução particular
xp_amp = gama / sqrt( (Wn^2 - We^2)^2 + (2*zeta*Wn*We)^2 );
phi = atan2( 2*zeta*Wn*We, (Wn^2 - We^2) );
xp = xp_amp * sin(We*T - phi);

% solução homogênea 
    if zeta == 0
    % Caso não amortecido
     C1 = y0(1) - xp(1);
     C2 = (y0(2) - We*xp_amp*cos(phi)) / Wn;
     xh = C1*cos(Wn*T) + C2*sin(Wn*T);

    elseif zeta < 1
    % Subamortecido
     Wd = Wn*sqrt(1 - zeta^2);
     C1 = y0(1) - xp(1);
     C2 = (y0(2) - We*xp_amp*cos(phi) + zeta*Wn*C1) / Wd;
     xh = exp(-zeta*Wn*T) .* ( C1*cos(Wd*T) + C2*sin(Wd*T) );

    elseif abs(zeta - 1) < 1e-6
    % Criticamente amortecido
     C1 = y0(1) - xp(1);
     C2 = y0(2) - We*xp_amp*cos(phi) + Wn*C1;
     xh = (C1 + C2*T) .* exp(-Wn*T);

    elseif zeta > 1
    % Superamortecido
     s1 = -Wn*(zeta - sqrt(zeta^2 - 1));
     s2 = -Wn*(zeta + sqrt(zeta^2 - 1));
     A = [1, 1; s1, s2] \ ([y0(1) - xp(1); y0(2) - We*xp_amp*cos(phi)]);
     C1 = A(1); C2 = A(2);
     xh = C1*exp(s1*T) + C2*exp(s2*T);
    end

x_t = (xh + xp).';     % para ser um vetor linha

% Erro absoluto

erro = abs(Y(1,:) - x_t);

% Gráficos


subplot(2,1,1);hold on; box on;
plot(T, Y(1,:), 'b', 'DisplayName','RK45 numérico');
plot(T, x_t, '--r', 'DisplayName','Analítico');
xlabel('t (s)'); ylabel('x(t)');
legend; title(['Oscilador amortecido, zeta = ', num2str(zeta)]);
grid on;

subplot(2,1,2);
plot(T, erro, 'k-','LineWidth',1);
xlabel('t (s)'); ylabel('Erro absoluto');
title('Erro da solução RK45 adaptativo');
grid on;

% Variação do tamanho do passo ao longo da simulação
figure; stairs(T(1:end-1), H, '-'); grid on;
xlabel('t (s) '); ylabel('dt'); title('Tamanho de passo adaptativo');