% Main Optimization File
clear all; clc;
x0 = [2.2 3.4 1.1];
A = [ ]; b = [ ];
Aeq = [ ]; beq = [ ];
lb = [ ]; ub = [ ];
nonlcon = [ ];

[x, fval] = fmincon(@objective_function,x0,A,b,Aeq,beq,lb,ub,nonlcon)