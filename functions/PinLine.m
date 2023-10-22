function ret = PinLine( x, y, z, ai, bi, ci, af, bf, cf )

%Essa função retorna "1" quando o ponto P = (x,y,z) pertence a linha que
%contém o segmento que começa em (ai,bi,ci) e termina em (af,bf,cf) e "0"
%caso contrário.

ret = 0;

A = [ai bi ci];
B = [af bf cf];
C = [x y z];

AB = A - B;
AC = A - C;

% if not collinear then return false as point cannot be on segment
if norm( cross( AB,AC ) ) < 0.0001
    ret = 1;
end

end