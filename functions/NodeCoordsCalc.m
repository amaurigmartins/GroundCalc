function [ node_coords ] = NodeCoordsCalc( coords_src )

% Essa função retorna uma matriz com as coordenadas de cada nó da malha de
% aterramento. As 3 primeiras colunas se referem as coordenadas x, y e z.
% Já a 4a coluna indica se o nó está no meio de algum segmento e qual é
% esse segmento. Quando o valor da 4a coluna é "0" isso indica que o nó não
% está no meio de nenhum segmento.

node_coords_rep = zeros((3*size(coords_src,1)),4);

j=1;

for i=1:size(coords_src,1)
    node_coords_rep(j,1) = coords_src(i,1);
    node_coords_rep(j,2) = coords_src(i,2);
    node_coords_rep(j,3) = coords_src(i,3);
    
    node_coords_rep(j+1,1) = coords_src(i,7);
    node_coords_rep(j+1,2) = coords_src(i,8);
    node_coords_rep(j+1,3) = coords_src(i,9);
    node_coords_rep(j+1,4) = i;                 %Ponto do meio do segmento i
    
    node_coords_rep(j+2,1) = coords_src(i,4);
    node_coords_rep(j+2,2) = coords_src(i,5);
    node_coords_rep(j+2,3) = coords_src(i,6);
    
    j = j + 3;
end

node_coords = unique(node_coords_rep,'rows', 'stable');

end

% Maneira de resolver sem usar a função unique
%
% check = ones(size(node_coords_rep,1),1);    %Para apagar linhas repetidas
% 
% for i=1:size(node_coords_rep,1)
%     for n=i+1:size(node_coords_rep,1)
%         if abs( node_coords_rep(i,1) - node_coords_rep(n,1) ) < tol
%             if abs( node_coords_rep(i,2) - node_coords_rep(n,2) ) < tol
%                 if abs( node_coords_rep(i,3) - node_coords_rep(n,3) ) < tol
%                     check(n) = 0;
%                 end
%             end
%         end
%     end
% end
% 
% node_coords = zeros(sum(check),4);
% 
% j = 1;
% 
% for i=1:size(node_coords_rep,1)
%     if check(i) ~= 0
%         node_coords(j,:) = node_coords_rep(i,:);
%         j = j+1;
%     end
% end
% 
% node_coords     %Matriz sem linhas repetidas