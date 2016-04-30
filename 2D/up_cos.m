function y = up_cos(t)


% a = [5e-1 5.538e-1 -4.621e-2 -8.663e-3 1.046e-3 -3.570e-4];

a = [0.500000000000000 ...
     0.553771275889086 ...
    -0.046202299134514 ...      
    -0.008656867928718 ...
     0.001044808241853 ...
    -0.000356152861644 ...
     0.000298368336856 ...
     0.000091404980607 ...
    -0.000006555455732 ...
    -0.000003383844607 ...
%      0.000012321122495
     ];

%% computation

omega = pi;
y = a(1) + a(2)*cos(omega*t) + a(3)*cos(3*omega*t) + ...
    a(4)*cos(5*omega*t) + a(5)*cos(7*omega*t) + a(6)*cos(9*omega*t) + ...
    a(7)*cos(11*omega*t) + a(8)*cos(13*omega*t) + a(9)*cos(15*omega*t); % + ...
%     a(11)*cos(17*omega*t);

y(abs(t)>=1) = 0;