clear all
global kalman_var;
kalman_var = zeros(0,8);
AA = diag([-12.4,-2.2,-2.2,-2.2]);

simout = sim('random');
simout = sim('test');