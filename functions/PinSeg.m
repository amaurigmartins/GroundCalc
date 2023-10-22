function ret = PinSeg( x, y, z, ai, bi, ci, af, bf, cf )

%Essa função retorna "1" quando o ponto P = (x,y,z) pertence ao segmento que
%começa em (ai,bi,ci) e termina em (af,bf,cf) e "0" caso contrário.

ret = 0;

A = [ai bi ci];
B = [af bf cf];
C = [x y z];

AB = A - B;
AC = A - C;

% if not collinear then return false as point cannot be on segment
if cross( AB,AC ) == 0
    % calculate the dotproduct of (AB, AC) and (AB, AB) to see point is now
    % on the segment
    dotAB = dot(AB, AB);
    dotAC = dot(AB, AC);
    % on end points of segment
    if dotAC == 0 || dotAC == dotAB
        ret = 1;
    % on segment
    elseif dotAC > 0 && dotAC < dotAB
        ret = 1;
    end
end

end