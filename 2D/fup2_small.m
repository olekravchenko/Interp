function y = fup2_small(x)


y = up_cos(0.25*x - 0.5) - 2*up_cos(0.25*x - 0.75) + 2*up_cos(0.25*x - 1) - ...
    2*up_cos(0.25*x - 1.25);

y = (36/13) * y;

y(abs(x)>2) = 0;
