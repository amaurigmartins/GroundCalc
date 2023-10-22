function [coords_src] = BuildMesh(width, length, depth, mesh_size)

% This function creates a grounding grid with the parameters chosen.

k = (width/mesh_size);
k2 = (length/mesh_size);

N_hor = k2+1;
N_ver = k+1; 

coords_src = zeros(( N_hor*k )+( N_ver*k2 ), 6);

for i=1:N_hor
   
    for j=1:k
        coords_src(j+(k*(i-1)),:) = [mesh_size*(j-1) (0 + (i-1)*mesh_size) -depth mesh_size*j (0 + (i-1)*mesh_size) -depth];
    end
    
end

for i=1:N_ver
    
    for j=1:k2
        coords_src(j+(k2*(i-1))+N_hor*k,:) = [(0 + (i-1)*mesh_size) mesh_size*(j-1) -depth (0 + (i-1)*mesh_size) mesh_size*j -depth];   
    end
    
end



end