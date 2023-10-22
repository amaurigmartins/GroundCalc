function [ Ii ] = CalcLeakageIi(coords_src, node_coords, V, Zg)

%This function calculates the currents dispersed in each segment of the
%grounding grid.

Ii = zeros(size(coords_src,1),1);
tol = 1E-5;

for i=1:size(coords_src,1)
    
    xm = coords_src(i,7);
    ym = coords_src(i,8);
    zm = coords_src(i,9);
    
    for n=1:size( node_coords, 1)
        
        if abs( xm - node_coords(n,2) )<tol
            if abs( ym - node_coords(n,3) )<tol
                if abs( zm - node_coords(n,4) )<tol
                    Node_m = node_coords(n,1);
                end
            end
        end
        
    end
    
    Ii(i) = V(Node_m)/Zg(i);
    
end

end