function xout = rk4(fun, dt, t0, x0)

f1 = fun(t0, x0);
f2 = fun(t0 + dt/2, x0 + (dt/2)*f1);
f3 = fun(t0 + dt/2, x0 + (dt/2)*f2);
f4 = fun(t0 + dt, x0 + dt*f3);

xout = x0 + (dt/6)*(f1+2*f2+2*f3+f4);