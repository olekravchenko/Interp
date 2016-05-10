function y = bspline( x )
% B-spline function of 3rd order

% x = x - 2;
y = zeros(size(x));
x0range = abs(x) < 1;
x1range = abs(x)>=1 & abs(x)<2;
x0 = x(x0range);
x1 = x(x1range);
c = 3/2;
y(x0range) = c * (2/3 - abs(x0).^2 + (1/2)*abs(x0).^3);
y(x1range) = c * (1/6)*(2-abs(x1)).*(2-abs(x1)).*(2-abs(x1));

% end