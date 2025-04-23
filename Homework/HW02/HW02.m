syms x1 x2 u1 u2 u3 s1 s2

eqn1 = 2*x1-2-2*u1-u2;
eqn2 = 2*x2-2-u1-2*u2;
eqn3 = -2*x1-x2+4+s1^2;
eqn4 = -x1 -2*x2+4+s2^2;

solution = solve([eqn1,eqn2,eqn3,eqn4],[x1, x2, u1, u2],"Real",true);
subs(solution.x1, [u1, s2],[0,0])
%%
Deqn1 = (-x1^2-x2^2-2*x1+16)*x1^2-18*x1*x2+13*x2^2-4;
Deqn2 = -x1^2-x2^2-2*x1+16;

H = hessian(Deqn1, [x1,x2]);

%% Problem 1 (a)
syms x1 x2 v1 

expr1 = -4*x1^2-3*x2^2+5*x1*x2+8;
eqConstr = x1+x2-4;
LM1 = expr1 + v1 * eqConstr;
LM_x1 = diff(LM1, x1)
LM_x2 = diff(LM1, x2)
LM_v1 = diff(LM1, v1)

sol1 = solve([LM_x1,LM_x2,LM_v1], [x1,x2,v1],'Real',true)

%% Problem 1 (b)
syms x1 x2 u1 v1 s1
expr2 = (x1 -2)^2+(x2 -1)^2;
ineqConstr2 = -x1-x2+4;
eqConstr2 = x1-x2-2;
LM2 = expr2+v1*eqConstr2+u1*(ineqConstr2+s1^2)
LM_x1 = diff(LM2, x1)
LM_x2 = diff(LM2, x2)
LM_v1 = diff(LM2, v1)
LM_u1 = diff(LM2, u1)

sol2 = solve([LM_x1,LM_x2,LM_v1,LM_u1], [x1,x2,v1,s1],"Real",true)

%% Problem 1 (c)

syms x1 x2 u1 s1
assume(u1 >= 0) % Lagrange multiplier for inequality constraint should be non-negative
assume(s1 >= 0) % 
expr3 = 9*x1^2-18*x1*x2+13*x2^2-4;
ineqConstr3 = -x1^2-x2^2-2*x1+16;

LM3 = expr3+u1*(ineqConstr3+s1^2)
LM_x1 = diff(LM3, x1)
LM_x2 = diff(LM3, x2)
LM_u1 = diff(LM3, u1)
LM_s1 = u1*s1
sol3 = vpasolve([LM_x1,LM_x2,LM_u1,LM_s1], [x1,x2,s1,u1])
%% Problem 2

A = [-1 3
    1 1
    1 -1
    -1 -3];
b = [10 6 2 -6];
f = [-2 -1];
[x,fval] = linprog(f,A,b)