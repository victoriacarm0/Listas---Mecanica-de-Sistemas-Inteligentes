clear
clc
close all

% Definição do sistema não suave (espaço de estados)
function dydt = sistema_naosuave(t,y,omega,zeta,wn,theta,A,Cp,R,g,zetab,wb)
    x  = y(1);  
    dx = y(2);  
    v  = y(3);
    
    if x < g
        ddx = -2*zeta*wn*dx - wn^2*x + theta*v + A*sin(omega*t);
    else
        % caso em que há contato com o batente
        ddx = -2*(zeta*wn + zetab*wb)*dx - wn^2*x - wb^2*(x - g) + theta*v + A*sin(omega*t);
    end
    
    dv = (-v/R - theta*dx)/Cp;
    
    dydt = [dx; ddx; dv];
end

% Função auxiliar: largura de banda (método 3dB)
function bw = calcBandwidth(freq, P)
    Pmax = max(P);
    if Pmax <= 0
        bw = 0; return;
    end
    P3dB = Pmax/2; % critério de meia potência
    idx = find(P >= P3dB); 
    if isempty(idx)
        bw = 0;
    else
        bw = freq(idx(end)) - freq(idx(1));
    end
end

% Parâmetros gerais
zeta = 0.025;         % amortecimento
wn   = 25;            % freq natural [rad/s]
theta= 0.0045;        % acoplamento N/V
Cp   = 4.2e-8;        % capacitância [F]
R    = 100e3;         % resistência [Ohm]
omega_range = 0:1:60; % variação de frequência excitação [rad/s]

% Parâmetros do modelo Biestável
alpha =  1;           
beta  = 1e4;          

% Parâmetros do modelo Não-suave (batente)
g_values = [0.001 0.005 0.01];  
zetab = 0.025;
wb = 5000;

% Diferentes amplitudes de excitação
A_values = [2.5 5 9.81];  

% Tempo de simulação
tspan = [0 20];        % segundos
transient = 3;        % tempo para desprezar transiente

% Função auxiliar para calcular potência média
calcPm = @(t,v) mean(v(t>transient).^2)/R;

% Loop para diferentes amplitudes A
for a_idx = 1:length(A_values)
    A = A_values(a_idx);

    % Inicializar resultados
    Pm_linear = zeros(size(omega_range));
    Pm_biestavel = zeros(size(omega_range));
    Pm_naosuave = zeros(length(g_values), length(omega_range));

    % Para salvar exemplos de v(t)
    v_examples = struct();

    % Diferentes frequências
    for k = 1:length(omega_range)
        omega = omega_range(k);

        % Sistema Linear
        f_linear = @(t,y)[ y(2);
            -2*zeta*wn*y(2) - wn^2*y(1) + theta*y(3) + A*sin(omega*t);
            ( -y(3)/R - theta*y(2) )/Cp ];
        [t,y] = ode45(f_linear, tspan, [0 0 0]);
        v_linear = y(:,3);
        Pm_linear(k) = calcPm(t,y(:,3));
        if k == round(length(omega_range)/2)
            v_examples.linear = struct('t',t,'v',y(:,3));
        end

        % Sistema Biestável 
        f_biestavel = @(t,y)[ y(2);
            -2*zeta*wn*y(2) + alpha*y(1) - beta*y(1)^3 + theta*y(3) + A*sin(omega*t);
            ( -y(3)/R - theta*y(2) )/Cp ];
        [t,y] = ode45(f_biestavel, tspan, [0 0 0]);
        v_biestavel = y(:,3);
        Pm_biestavel(k) = calcPm(t,y(:,3));
        if k == round(length(omega_range)/2)
            v_examples.biestavel = struct('t',t,'v',y(:,3));
        end

        % Sistema Não-suave (batente)
        for g_idx = 1:length(g_values)
            g = g_values(g_idx);
            f_naosuave = @(t,y) sistema_naosuave(t,y,omega,zeta,wn,theta,A,Cp,R,g,zetab,wb);
            [t,y] = ode45(f_naosuave, tspan, [0 0 0]);
            v_nsuave = y(:,3);
            Pm_naosuave(g_idx,k) = calcPm(t,y(:,3));
            if k == round(length(omega_range)/2)
                v_examples.naosuave{g_idx} = struct('t',t,'v',y(:,3));
            end
        end
    end



    % Cálculo da largura de banda 
    BW_linear = calcBandwidth(omega_range,Pm_linear);
    BW_biestavel = calcBandwidth(omega_range,Pm_biestavel);
    BW_naosuave = arrayfun(@(i) calcBandwidth(omega_range,Pm_naosuave(i,:)),1:length(g_values));

    % Mostrar resultados na tela
    fprintf('\nAmplitude A = %.2f\n',A);
    fprintf('  Linear     -> BW = %.2f rad/s\n',BW_linear);
    fprintf('  Biestavel  -> BW = %.2f rad/s\n',BW_biestavel);
    for g_idx = 1:length(g_values)
        fprintf('  Não-suave (g=%.3f) -> BW = %.2f rad/s\n', g_values(g_idx), BW_naosuave(g_idx));
    end

    % Plotagem dos resultados para essa amplitude A
    figure; hold on; grid on;
    p1 = plot(omega_range,Pm_linear,'LineWidth',1.5);
    p2 = plot(omega_range,Pm_biestavel,'LineWidth',1.5);
    p_nao = gobjects(length(g_values),1);
    for g_idx = 1:length(g_values)
        p_nao(g_idx) = plot(omega_range,Pm_naosuave(g_idx,:),'LineWidth',1.5);
    end
    xlabel('\omega [rad/s]');
    ylabel('P_m [W]');
    title(['Potência média x Frequência - A = ' num2str(A)]);
    leg = [{'Linear','Biestável'}, ...
           arrayfun(@(g) sprintf('Não-suave (g=%.3f)',g), g_values, 'UniformOutput',false)];
    legend([p1 p2 p_nao(:)'], leg, 'Location','NorthEast');

    % Plotar os sinais de v(t) para comparação
    figure; hold on; grid on;
    plot(v_examples.linear.t,v_examples.linear.v,'LineWidth',1.5);
    plot(v_examples.biestavel.t,v_examples.biestavel.v,'LineWidth',1.5);
    for g_idx = 1:length(g_values)
        plot(v_examples.naosuave{g_idx}.t,v_examples.naosuave{g_idx}.v,'LineWidth',1.5);
    end
    xlabel('Tempo [s]');
    ylabel('v(t) [V]');
    title(['Resposta em tensão v(t) - A = ' num2str(A)]);
    legend(leg,'Location','Best');

end
