function [Zint] = CalcZint(coords_src,f,cond_rho,cond_mur,cond_radius)

%This function calculates the internal impedance of each segment in the
%grounding grid.

% rho_cond = 1.72E-8; %Resistividade elétrica do cobre
u_0 = 4*pi*1e-7; %Permeabilidade magnética do vácuo

Zint = zeros(size(coords_src,1),1);
l = coords_src(:,10); %Cria um vetor com os comprimentos de cada segmento
% cond_radius = coords_src(:,11); %Cria um vetor com os raios de cada segmento

for n=1:size(coords_src,1)

    R = ( cond_rho(n)*l(n) )/( pi*(cond_radius(n)^2) ); %Fórmula (2.14) do pdf "Transient-behavior..."
    L = ( (u_0*cond_mur(n))/(2*pi) )*( log( 2*l(n)/cond_radius(n) ) - 1 ); %Fórmula (2.15) do pdf "Transient-behavior..."

    Zint(n) = complex(R,2*pi*f*L); %Impedância interna do segmento n

end

% Fórmula do pdf Zint

% u_cond=1;
% 
% sigma_i = sqrt( rho_cond/pi*60*u_cond );
% 
% q = (2*cond_radius(1))/(sqrt(2)*sigma_i);
% 
% Rs = sqrt(pi*60*u_cond*rho_cond);
% 
% Zint = complex( 0 , (l(1)*Rs)/(sqrt(2)*pi*cond_radius(1)) )*( (complex(be)) )


end