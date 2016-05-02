function [M, rhs] = AssembleMatrix(method, nx, hx, yai, d1yai, gv, d1gv)

ind = 2:nx-1;       % Index array
npt = nx-2;         % Number of inner points

switch method
    case 'IDO'
        % IDO scheme
        % matrix
        M11             =    (2/hx^2)*(diag(ones(1,nx-3),-1)+diag(ones(1,nx-3),1)-2*diag(ones(1,nx-2)));
        M12             =    (0.5/hx)*(diag(ones(1,nx-3),-1)-diag(ones(1,nx-3),1));
        M21             = -(7.5/hx^3)*(diag(ones(1,nx-3),-1)-diag(ones(1,nx-3),1));
        M22             = -(1.5/hx^2)*(diag(ones(1,nx-3),-1)+diag(ones(1,nx-3),1)+8*diag(ones(1,nx-2)));
        M               = vertcat(horzcat(M11,M12), horzcat(M21,M22));
        % right hand side
        rhs             = [gv(ind) d1gv(ind)];
        rhs(1)          = rhs(1)        - (2/hx^2)*yai(1)       - (0.5/hx)*d1yai(1);
        rhs(nx-2)       = rhs(nx-2)     - (2/hx^2)*yai(2)       + (0.5/hx)*d1yai(2);
        rhs(nx-1)       = rhs(nx-1)     + (7.5/hx^3)*yai(1)     + (1.5/hx^2)*d1yai(1);
        rhs(2*(nx-2))   = rhs(2*(nx-2)) - (7.5/hx^3)*yai(2)     + (1.5/hx^2)*d1yai(2);
    case 'CIPBS0'
        % CIPBS0 scheme
        % matrix
        M               = (1/hx)*diag(ones(1,nx-3),-1) -(2/hx)*eye(nx-2) + (1/hx)*diag(ones(1,nx-3),1);
        % right hand side vector (2:nx-1)*1 size        
%         rhs             = zeros(nx-2,1);
        rhs             = (hx/6)*gv(ind-1)+(2*hx/3)*gv(ind)+(hx/6)*gv(ind+1);
        rhs(1)          = rhs(1)        - yai(1)/hx;
        rhs(nx-2)       = rhs(nx-2)     - yai(2)/hx;
    case 'CIPBS1'
        % CIPBS1 scheme
        M1 = zeros(npt, 2*npt);
        M2 = M1;
        % matrix M (1:nx-2)*(1:nx-2) size
        for k=2:npt-1
            M1(k,2*k-3:2*k+2)   = [6/(5*hx) 1/10 -12/(5*hx) 0 6/(5*hx) -1/10];
            M2(k,2*k-3:2*k+2)   = [-1/10 hx/30 0 -(4*hx)/15 1/10 hx/30];
        end
        M1(1,1:4)               = [-12/(5*hx) 0 6/(5*hx) -1/10];
        M1(npt,2*npt-3:2*npt)   = [6/(5*hx) 1/10 -12/(5*hx) 0];
        M2(1,1:4)               = [0 -(4*hx)/15 1/10 hx/30];
        M2(npt,2*npt-3:2*npt)   = [-1/10 hx/30 0 -(4*hx)/15];
        
        M = vertcat(M1, M2);
        
        % right hand side vector (1:nx-2)*1 size
        rhs1        = (9*hx/70)*gv(ind-1)+(26*hx/35)*gv(ind)+(9*hx/70)*gv(ind+1) + ...
                      (13*hx^2/420)*d1gv(ind-1)-(13*hx^2/420)*d1gv(ind+1);
        rhs1(1)     = rhs1(1)    - (6/(5*hx)*yai(1) + (1/10)*d1yai(1));
        rhs1(nx-2)  = rhs1(nx-2) - (6/(5*hx)*yai(2) - (1/10)*d1yai(2));
        % second part of right hand side vector (2:nx-1)*1 size
        rhs2        = -(13*hx^2/420)*gv(ind-1)+(13*hx^2/420)*gv(ind+1) - ...
                       (hx^3/140)*d1gv(ind-1)+(2*hx^3/105)*d1gv(ind)-(hx^3/140)*d1gv(ind+1);
        rhs2(1)     = rhs2(1)    - ((-1/10)*yai(1) + (hx/30)*d1yai(1));
        rhs2(nx-2)  = rhs2(nx-2) - ((1/10)*yai(2) + (hx/30)*d1yai(2));
        
        rhs = [rhs1 rhs2];
    otherwise
        disp('Warning: Incorrect Method');
        M = 0; rhs = 0;
        return
end % switch