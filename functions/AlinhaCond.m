function [coords_src] = AlinhaCond(coords_src)

%This function makes all vertical conductors be directed down.

tol = 1E-5;

for i=1:size(coords_src,1)
    
    zs = coords_src(i,3);
    ze = coords_src(i,6);
    
    if abs(zs-ze)>tol % Se o segmento é vertical
        
        if zs<ze
           
            z_aux = zs;
            zs = ze;
            ze = z_aux;
            
        end
        
    end
    
    coords_src(i,3) = zs;
    coords_src(i,6) = ze;
    
end

end