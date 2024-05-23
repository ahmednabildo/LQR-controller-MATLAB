clear
clc
close all

%Define system

A = [0 1 0 0;0 0 1 0;0 0 0 1;0 0 0 -5];

B = [0;0;0;17];

C = [1 0 0 0];

D = 0;

%Choose Q and R

Q = diag([1 0.5 0.2 0.1]);
R = [0.000005];
        

[K, S, E] = lqr(A, B, Q, R);

disp('K computed via LQR:') 
K
Acl = (A - B*K);
Ecl = eig(Acl);
syscl = ss(Acl,B,C,D);
Kdc = dcgain(syscl);
Kr = 1 / Kdc;
syscl_sc = ss(Acl,B*Kr,C,D);
[y, t] = step(syscl_sc);
step(syscl_sc);
[num , den] = ss2tf(Acl,B*Kr,C,D);
G = tf(num , den);

[num1 , den1] = ss2tf(A,B,C,D);
P = tf(num1 , den1);
C = G/(P-G*P);
% Calculate step response characteristics using stepinfo
info = stepinfo(y,t);
disp('info computed via LQR:') 
info

