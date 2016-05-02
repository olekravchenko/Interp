function [out0, out1, out2] = InterpPoissonSolver(method, f, varargin)

if length(varargin) == 0
    a   = 0;
    b   = 1;
    nx  = 10;
elseif length(varargin) == 2
    a   = varargin{1};
    b   = varargin{2};
    nx  = 10;
elseif  length(varargin) == 3
    a   = varargin{1};
    b   = varargin{2};
    nx  = varargin{3};
else
    disp('Incorrect Input');
    out0 = 0; out1 = 0; out2 = 0;   
    return
end % if

% Space grid
out0 = linspace(a,b,nx);
hx   = out0(2) - out0(1);

%% Symbolic functions
tic
yafunc      = sym(f);
d1yafunc    = diff(yafunc);
d2yafunc    = diff(d1yafunc);
gfunc       = d2yafunc;
d1gfunc     = diff(gfunc);
g           = matlabFunction(gfunc);
d1ya        = matlabFunction(d1yafunc);
d1g         = matlabFunction(d1gfunc);
t = toc;

disp(['Symbolic Operations: ' num2str(t)])

%% Boundary conditions
% gi      = [g(out0(1))       g(out0(nx))];
yai     = [f(out0(1))       f(out0(nx))];
% d1gi    = [d1g(out0(1))     d1g(out0(nx))];
d1yai   = [d1ya(out0(1))    d1ya(out0(nx))];

%% Vectors
gv      = g(out0);
d1gv    = d1g(out0);

switch method
    case 'IDO'
        [M, rhs] = AssembleMatrix(method,nx,hx,yai,d1yai,gv,d1gv);       
        % Solve linear system
        coeff = M \ rhs';
        % coeff = bicgstab(M, rhs');
    case 'CIPBS0'
        [M, rhs] = AssembleMatrix(method,nx,hx,yai,d1yai,gv,d1gv);
        % Solve linear system
        coeff = M \ rhs';
        % Browse coeff vector
        coeff = vertcat(coeff,zeros(2*(nx-2),1));
    case 'CIPBS1'
        [M, rhs] = AssembleMatrix(method,nx,hx,yai,d1yai,gv,d1gv);
        % Solve linear system
        coeff = M \ rhs';
        coeff = [coeff(1:2:end)' coeff(2:2:end)']';
        % coeff = bicgstab(M, rhs');
    otherwise
        disp('Warning: Incorrect Method');
        out0 = 0; out1 = 0; out2 = 0;
        return
end % switch

% Output
out1 = [yai(1) coeff(1:nx-2)' yai(2)];
out2 = [d1yai(1) coeff(nx-1:2*(nx-2))' d1yai(2)];