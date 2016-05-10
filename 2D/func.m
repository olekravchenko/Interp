function [output, info] = func(id,x,y)

switch id
    case '0'
        output  = 1 + 0 * (x .* y);        
        info    = '$$1$$';        
    case '1'
        output  = 0.05 * (x + y);
        info    = '$$(x + y) / 200$$';
    case '2'
        output  = 0.01 * x .* y;
        info    = '$$xy/100$$';
    case '2a'
        output  =  0.0069 * (x - 3).^2 - y.^2 *0.01;          
        info    = '$$(x-3)^2/144 - y^2/100$$';                       
    case '3'
        output  = (x.^2 - x.*y) / 200;
        info    = '$$(x^2 - xy) / 200$$';        
    case '4'
        output  = x.*(10-x).*y.*(10-y) / 625;        
        info    = '$$xy(10-x)(10-y) / 625$$';
    case '5'        
        output  = x.^2 .* y.^2 * 1e-4;        
        info    = '$$x^2y^2 10^{-4}$$';
    case '6'
        output  = sin(0.2 * x .* y);                
        info    = '$$\sin(0.2xy)$$';
    case '7'
        output  = sin(0.2 * pi * x) .* sin(0.25 * pi * y);
        info    = '$$\sin(0.2\pi x)\sin(0.25\pi y)$$';
    case '8'
        output  = cos(0.25 * pi * x) .* cos(0.2 * pi * y);                        
        info    = '$$\cos(0.25 \pi x)\cos(0.2\pi y)$$';
    case '9'
        output  = (x.^2 + y.^2 - 1) .* sin(0.2*x.*y) / 200;          
        info    = '$$(x^2 + y^2 - 1)\sin(0.2xy)/200$$';
    case '10'
        output  = 1 + tanh(10*x.^2.*y^2);          
        info    = '$$1 + \tanh(10x^2y^2)$$';        
    case '11'
        output  = x.*y.^3 * 1e-5;          
        info    = '$$xy^3 10^{-5}$$';
    case '12'
        output  = x.^2+y.^2;          
        info    = '$$x^2 + y^2$$';
    case '13' 
        output  = -2.*sin(x).*sin(y).+0*(x .* y);
        info    = '-2sin(x)*sin(y)';
    case '14' 
        output  =  -8.*sin(2.*x).*sin(2.*y);
        info    = '-2sin(x)*sin(y)';
    case '15' 
        output  =  -16.*pi*sin(2.*pi*x).*sin(2.*pi*y);
        info    = '-2sin(x)*sin(y)';
    case '16' 
        output  =  2.*(x.^2+y.^2-x-y);
        info    = '-2sin(x)*sin(y)';
    otherwise
        output  = 1e-5 + 0*(x .* y);
        info    = '$$10^{-5}$$';
end