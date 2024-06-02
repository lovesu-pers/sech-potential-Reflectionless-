function [psi, E, conv] = egenfunc(Nx, Ux, m, hbar, dx, n)
% 
%     h1 = ones(Nx,1)*(hbar^2)/(2*m*dx^2) + Ux;
%     h2 = -ones(Nx,1)*(hbar^2)/(m*dx^2);
%     
%     H = spdiags([h2 h1 h2], -1:1, Nx, Nx);
%     
      C = -hbar^2/(2*m*dx^2);
    
      h1 = ones(Nx,1) .* (-C*30/(12) + Ux);
      h2 = ones(Nx,1) * (C*16/(12));
      h3 = ones(Nx,1) * (-C/(12));
      H = spdiags([h3 h2 h1 h2 h3], -2:2, Nx, Nx);
      
      H(end,1) = -hbar/(2*m*dx^2);
      H(1,end) = -hbar/(2*m*dx^2);
%       H(1,1:3) = C*[1, -2, 1];
%       H(1,(Nx-4):(Nx-2)) = C*[1, -2, 1];

      [psi, E, conv] = eigs(H, n, 'smallestreal');
end