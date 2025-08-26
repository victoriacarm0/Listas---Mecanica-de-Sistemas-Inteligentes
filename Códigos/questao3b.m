clear 
close all
clc

% Pêndulo simples não amortecido
% Definição simbólica
syms theta_1 theta_2 Wn real

% Vetor de estados
x = [theta_1; theta_2];

% Sistema em espaço de estados
f = [ theta_2;
     -Wn^2*sin(theta_1)];

% Determinar os pontos de equilíbrio
eqs = f == [0; 0];
sol_eq = solve(eqs, [theta_1, theta_2], 'Real', true);
equilibrio = [sol_eq.theta_1, sol_eq.theta_2];
disp('Pontos de equilíbrio');
disp(equilibrio);

% Calcular a Jacobiano 
J = jacobian(f, x);
disp('Jacobiano');
disp(J);

% Autovalores gerais da Jacobiana
syms lambda
autovalores_gerais = eig(J);

disp('Autovalores gerais');
autovalores_gerais = simplify(autovalores_gerais);
disp(autovalores_gerais);

% Armazenar os autovalores de cada ponto em uma célula
autovalores_por_ponto = cell(length(equilibrio), 1);

% Avaliar os autovalores em cada ponto de equilíbrio e armazenar
 for k = 1:length(equilibrio)
    % substitui ponto de equilíbrio nos autovalores
    autovalores_k = subs(autovalores_gerais, ...
        [theta_1, theta_2], [equilibrio(k,1), equilibrio(k,2)]);
    
    % salva na célula
    autovalores_por_ponto{k} = simplify(autovalores_k);
    
    % mostra resultado simbólico
    fprintf('\nPonto de equilíbrio %d: (%s, %s)\n', ...
        k, char(equilibrio(k,1)), char(equilibrio(k,2)));
    disp('Autovalores:')
    disp(autovalores_por_ponto{k})
 end

% Avaliação da natureza de estabilidade dos pontos
% Valores numéricos de exemplo
Wn = 2*pi;

for k = 1:length(autovalores_por_ponto)
    % substitui valores numéricos
    aut_num = double(subs(autovalores_por_ponto{k}, ...
        [alpha, beta, zeta], [alpha_val, beta_val, zeta_val]));
    
    fprintf('\n>> Ponto de equilíbrio %d\n', k);
    disp('Autovalores numéricos:')
    disp(aut_num)
    
    % Avaliação 
    if all(real(aut_num) < 0) && all(imag(aut_num) == 0)
        disp('Ponto tipo sorvedouro (estável)');
    elseif all(real(aut_num) > 0) && all(imag(aut_num) == 0)
        disp('Ponto tipo fonte (instável)');
    elseif (real(aut_num(1)) * real(aut_num(2)) < 0)
        disp('Ponto tipo sela (instável)');
    elseif all(real(aut_num) == 0) && any(imag(aut_num) ~= 0)
        disp('Ponto tipo centro (neutro)');
    elseif all(real(aut_num) < 0) && any(imag(aut_num) ~= 0)
        disp('Ponto espiral estável (foco)');
    elseif all(real(aut_num) > 0) && any(imag(aut_num) ~= 0)
        disp('Ponto espiral instável (foco)');
    else
        disp('Outro tipo de ponto ou caso não tratado.');
    end
end