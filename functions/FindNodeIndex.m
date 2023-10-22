function [node_index] = FindNodeIndex(x_node, y_node, z_node, node_coords)

% This function returns the node index of given coordinates

node_index = 0;

for i=1:size(node_coords)
    
    x = node_coords(i,2);
    y = node_coords(i,3);
    z = node_coords(i,4);
    
    if (x==x_node) && (y==y_node) && (z==z_node)
        node_index = i;
    end
    
end

if node_index == 0
    fprintf("There is no node in the coordinates given.")
end

end