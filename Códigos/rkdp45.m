function [yout4, yout5] = rkdp45(fun, dt, t0, y0)

    % Coeficientes do Dormand–Prince 4(5)
    c2 = 1/5;
    c3 = 3/10;
    c4 = 4/5;
    c5 = 8/9;
    c6 = 1;
    c7 = 1;

    a21 = 1/5;
    a31 = 3/40;        a32 = 9/40;
    a41 = 44/45;       a42 = -56/15;      a43 = 32/9;
    a51 = 19372/6561;  a52 = -25360/2187; a53 = 64448/6561;  a54 = -212/729;
    a61 = 9017/3168;   a62 = -355/33;     a63 = 46732/5247;  a64 = 49/176;     a65 = -5103/18656;
    a71 = 35/384;      a72 = 0;           a73 = 500/1113;    a74 = 125/192;    a75 = -2187/6784;   a76 = 11/84;

    % Pesos da 5ª ordem
    b1 = 35/384; b2 = 0; b3 = 500/1113; b4 = 125/192; b5 = -2187/6784; b6 = 11/84; b7 = 0;

    % Pesos da 4ª ordem
    b1s = 5179/57600; b2s = 0; b3s = 7571/16695; b4s = 393/640;
    b5s = -92097/339200; b6s = 187/2100; b7s = 1/40;

    % Cálculo das etapas
    k1 = fun(t0,               y0);
    k2 = fun(t0 + c2*dt,       y0 + dt*(a21*k1));
    k3 = fun(t0 + c3*dt,       y0 + dt*(a31*k1 + a32*k2));
    k4 = fun(t0 + c4*dt,       y0 + dt*(a41*k1 + a42*k2 + a43*k3));
    k5 = fun(t0 + c5*dt,       y0 + dt*(a51*k1 + a52*k2 + a53*k3 + a54*k4));
    k6 = fun(t0 + c6*dt,       y0 + dt*(a61*k1 + a62*k2 + a63*k3 + a64*k4 + a65*k5));
    k7 = fun(t0 + c7*dt,       y0 + dt*(a71*k1 + a72*k2 + a73*k3 + a74*k4 + a75*k5 + a76*k6));

    % Estimativa de 5ª ordem
    yout5 = y0 + dt*(b1*k1 + b2*k2 + b3*k3 + b4*k4 + b5*k5 + b6*k6 + b7*k7);

    % Estimativa de 4ª ordem
    yout4 = y0 + dt*(b1s*k1 + b2s*k2 + b3s*k3 + b4s*k4 + b5s*k5 + b6s*k6 + b7s*k7);
end