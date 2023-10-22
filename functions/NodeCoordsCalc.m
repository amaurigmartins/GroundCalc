function [ node_coords ] = NodeCoordsCalc( coords_src )

% Essa função retorna uma matriz com as coordenadas de cada nó da malha de
% aterramento. As colunas se referem as coordenadas x, y e z.

node_coords_rep = zeros((3*size(coords_src,1)),4);

j=1;

for i=1:size(coords_src,1)
    node_coords_rep(j,2) = coords_src(i,1);
    node_coords_rep(j,3) = coords_src(i,2);
    node_coords_rep(j,4) = coords_src(i,3);
    
    node_coords_rep(j+1,2) = coords_src(i,7);
    node_coords_rep(j+1,3) = coords_src(i,8);
    node_coords_rep(j+1,4) = coords_src(i,9);
    
    node_coords_rep(j+2,2) = coords_src(i,4);
    node_coords_rep(j+2,3) = coords_src(i,5);
    node_coords_rep(j+2,4) = coords_src(i,6);
    
    j = j + 3;
end

node_coords = unique(node_coords_rep,'rows', 'stable');

% node_coords(size(node_coords,1),:) = []; # Apaga a última linha

for i=1:size(node_coords,1)
    node_coords(i,1) = i;
end

end