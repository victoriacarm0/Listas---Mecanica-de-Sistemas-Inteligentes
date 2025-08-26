clear 
close all
clc

% Sistema multiestável com 2 GDL
% Definição simbólica
syms x1 x2 v1 v2 real

% Valores das constantes (essa solução será inteiramente numérica)
z1 = 0.05;
z2 = 0.08;
a1 = -2;
a2 = -2;
b1 = 1;
b2 = 1.5;
rho = 0.5;
Os = 3 ;

% Vetor de estados
X = [x1; v1; x2; v2];

% Sistema em espaço de estados
f = [ v1;
      -2*z1*v1 + 2*z2*(v2 - v1) - (1+a1)*x1 - b1*x1^3 + rho*Os^2*(x2-x1);
      v2;
     (-(2*z2)*(v2 - v1) - a2*x2 - b2*x2^3 - rho*Os^2*(x2 - x1))/rho];

% Determinar os pontos de equilíbrio
eqs = f == 0;
sol_eq = solve(eqs, X, 'Real', true);
equilibrio = double([sol_eq.x1, sol_eq.v1, sol_eq.x2, sol_eq.v2]);
disp('Pontos de equilíbrio');
disp(equilibrio);

% Calcular a Jacobiano 
J = jacobian(f, X);
disp('Jacobiano');
disp(J);

% Autovalores gerais da Jacobiana
 [V, autovalores_gerais] = eig(J);
 
% disp('Autovalores gerais');
autovalores_gerais = simplify(autovalores_gerais);

% Armazenar os autovalores de cada ponto em uma célula
autovalores_por_ponto = cell(length(equilibrio), 1);

% Avaliar os autovalores em cada ponto de equilíbrio e armazenar
 for k = 1:size(equilibrio, 1)
    % substitui ponto de equilíbrio nos autovalores
    autovalores_k = subs(autovalores_gerais, ...
        [x1, v1, x2, v2], [equilibrio(k,1), equilibrio(k,2), equilibrio(k,3), equilibrio(k,4)]);
    
    % salva na célula
    autovalores_por_ponto{k} = simplify(autovalores_k);
    
 end

% Avaliação da natureza de estabilidade dos pontos
for k = 1:length(autovalores_por_ponto)
    aut_sym = autovalores_por_ponto{k};     
    aut_num = diag(double(aut_sym));       % muda de simbólico para numérico

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