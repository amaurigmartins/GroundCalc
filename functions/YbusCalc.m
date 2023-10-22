function [ Ybus ] = YbusCalc( node_coords, coords_src, Y_int, Yg, tol, l )

%Essa função retorna a matriz de admitâncias Ybus do crcuito equivalente da malha de aterramento.

%%% Cálculo das admitâncias Y_conex referentes as conexões externas de cada nó

Y_conex = zeros(size(node_coords,1),1);

for i=1:size(node_coords,1)
    
    this_x = node_coords(i,1);
    this_y = node_coords(i,2);
    this_z = node_coords(i,3);
    
    for n=1:size(coords_src,1)

        if abs( this_x-coords_src(n,1) ) < tol                      %O segmento n começa no nó i
            if abs( this_y-coords_src(n,2) ) < tol
                if abs( this_z-coords_src(n,3) ) < tol
                    Y_conex(i) = Y_conex(i) + (Y_int(n)/2);
                    fprintf('O nó %i está no começo do segmento %i!\n',i,n);
                end
            end
        end
        if abs( this_x-coords_src(n,4) ) < tol                  %O segmento n termina no nó i
            if abs( this_y-coords_src(n,5) ) < tol
                if abs( this_z-coords_src(n,6) ) < tol
                    Y_conex(i) = Y_conex(i) + (Y_int(n)/2);
                    fprintf('O nó %i está no final do segmento %i!\n',i,n);
                end
            end
        end
        if abs( this_x-coords_src(n,7) ) < tol                  %O nó i está no meio do segmento n
            if abs( this_y-coords_src(n,8) ) < tol
                if abs( this_z-coords_src(n,9) ) < tol
                    Y_conex(i) = Y_conex(i) + Y_int(n) + Yg(n);
                    fprintf('O nó %i está no meio do segmento %i!\n',i,n);
                end
            end
        end
    end
    
end

% Y_conex

%%% Construindo a matriz de admitâncias Ybus

Ybus = zeros(size(node_coords,1));

for i=1:(size(node_coords,1))
    
    for j=i:(size(node_coords,1))
        
        if i==j %Na diagonal principal
            
            Ybus(i,i) = Y_conex(i); 
        
        else   %Fora da diagonal principal
            
            xi = node_coords(i,1);
            yi = node_coords(i,2);
            zi = node_coords(i,3);
            
            xj = node_coords(j,1);
            yj = node_coords(j,2);
            zj = node_coords(j,3);
            
            dist_ij = sqrt( (xj-xi)^2 + (yj-yi)^2 + (zj-zi)^2 );

            if node_coords(i,4)     %O nó i está no meio do segmento n
                n = node_coords(i,4);
                if abs( dist_ij - (l(n)/2) ) < tol  %E a distancia entre os nós é metade do segmento n
                    Ybus(i,j) = -Y_int(n)/2;
                end
            elseif node_coords(j,4)     %O nó j está no meio do segmento n
                n = node_coords(j,4);
                if abs( dist_ij - (l(n)/2) ) < tol  %E a distancia entre os nós é metade do segmento n
                    Ybus(i,j) = -Y_int(n)/2;
                end
            end
            
        end

    end
    
end

for i=1:size(Ybus,1)    %Completa o que está abaixo da diagonal pois a matriz é simétrica
    for j=(i+1):size(Ybus,1)
        Ybus(j,i) = Ybus(i,j);
    end
end

Ybus

end