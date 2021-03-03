function [ soma ] = greenfun2la( rho1, rho2, h, xs, ys, zs, xe, ye, ze, comp, x0, y0, z0 )

TOL=0.0001;
k=(rho2-rho1)/(rho2+rho1);
e=abs(zs);
zs=abs(zs);
ze=abs(ze);
z0=abs(z0);
% l=sqrt((xs-xe)^2+(ys-ye)^2+(zs-ze)^2);
l=comp;
if abs(zs-ze) <= TOL % horizontal
    alfa=atan2((ye-ys),(xe-xs));
    A=[cos(alfa) -sin(alfa); sin(alfa) -cos(alfa)];
    B=[x0-xs; y0-ys];
    rot=[cos(alfa) sin(alfa); -sin(alfa) cos(alfa)];
    transf=rot*B;
%     if rcond(A)<1e-6
%         transf=pinv(A)*B;
%     else
%         transf=linsolve(A,B);
%     end
    u0=transf(1,1);
    v0=transf(2,1);
    w0=z0-zs;
else % vertical
    u0=z0-zs;
    v0=ys-y0;
    w0=xs-x0;
end
% gama=@(x) log(((l-u0)+sqrt((l-u0)^2+v0^2+(2*x+w0)^2))/(-u0+sqrt(u0^2+v0^2+(2*x+w0)^2)));
psi=@(x) log((l - u0 + (v0^2 + w0^2 + (2*x + w0)^2 + (l - u0)^2)^(1/2))/((u0^2 + v0^2 + w0^2 + (2*x + w0)^2)^(1/2) - u0));
s=psi(TOL)+psi(e);
M = 10;
n = 1;

if zs <= h % src at top layer
    if z0 <= h % obs at top layer
        %Equação 5-2
        A=(rho1)/(4*pi);
        while abs(M) > TOL
            M = k^n*(psi(n*h)+psi(n*h+e)+psi(-n*h)+psi(-n*h+e));
            s = s + M;
            n = n + 1;
        end
    else % ob at bottom layer
        %eqn 5-3
        A=(rho1*(1+k))/(4*pi);
        while abs(M) > TOL
            M = k^n*(psi(n*h)+psi(n*h+e));
            s = s + M;
            n = n + 1;
        end
    end
else % src at bottom layer
    if z0 <= h %obs at top layer
        %eqn 5-4
        A=(rho2)/(4*pi);
        while abs(M) > TOL
            M = k^n*(psi(-n*h)+psi(n*h+e)-psi((-n-1)*h)-psi((n-1)*h+e));
            s = s + M;
            n = n + 1;
        end
    else %obs at bottom layer
        %eqn 5-5
        A=(rho2)/(4*pi);
        while abs(M) > TOL
            M = k^n*(psi(n*h+e)-psi((n-2)*h+e));
            s = s + M;
            n = n + 1;
        end
    end
end
soma = A*s;
% if isnan(soma) || isinf(soma)
%     keyboard
% end
end

