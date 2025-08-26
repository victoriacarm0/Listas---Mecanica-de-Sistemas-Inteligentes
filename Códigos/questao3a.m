 clear 
close all
clc

% Oscilador de Duffing 
% Definição simbólica
syms x_1 x_2 alpha beta zeta real

% Vetor de estados
x = [x_1; x_2];

% Sistema em espaço de estados
f = [ x_2;
     -2*zeta*x_2 - alpha*x_1 - beta*x_1^3 ];

% Determinar os pontos de equilíbrio
eqs = f == [0; 0];
sol_eq = solve(eqs, [x_1, x_2], 'Real', true);
equilibrio = [sol_eq.x_1, sol_eq.x_2];
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
        [x_1, x_2], [equilibrio(k,1), equilibrio(k,2)]);
    
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
alpha_val = -1.2; beta_val = 0.3; zeta_val = 0.05;

for k = 1:length(autovalores_por_ponto)
    % substitui valores numéricos
    aut_num = double(subs(autovalores_por_ponto{k}, ...
        [alpha, beta, zeta], [alpha_val, beta_val, zeta_val]));
    
    fprintf('\n>> Ponto de equilíbrio %d\n', k);
    disp('Autovalores numéricos:')
    disp(aut_num)
    
    % Avaliação 
    if all(real(aut_num) < 0) && all(imag(aut_num) == 0)
        disp('Ponto tipo sorvedouro (nó estável)');
    elseif all(real(aut_num) > 0) && all(imag(aut_num) == 0)
        disp('Ponto tipo fonte (instável)');
    elseif (any(real(aut_num) > 0) && any(real(aut_num) < 0))
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