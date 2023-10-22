function [NetList] = CalcNetList(coords_src, node_coords, Zg, Zint)

%This function returns the NetList of the grounding system equivalent
%circuit. It considers the "T" equivalent model for each segment.

NetList = zeros( 3*size(coords_src,1), 4 );

tol = 1E-5;
j = 1;

for i=1:size( coords_src,1 )
    
    xi = coords_src(i,1);
    yi = coords_src(i,2);
    zi = coords_src(i,3);
    
    xf = coords_src(i,4);
    yf = coords_src(i,5);
    zf = coords_src(i,6);
    
    xm = coords_src(i,7);
    ym = coords_src(i,8);
    zm = coords_src(i,9);
    
    for n=1:size( node_coords,1 )
       
        if abs( xi - node_coords(n,2) )<tol
            if abs( yi - node_coords(n,3) )<tol
                if abs( zi - node_coords(n,4) )<tol
                    Node_i = node_coords(n,1);
                end
            end
        end
        
        if abs( xm - node_coords(n,2) )<tol
            if abs( ym - node_coords(n,3) )<tol
                if abs( zm - node_coords(n,4) )<tol
                    Node_m = node_coords(n,1);
                end
            end
        end
        
        if abs( xf - node_coords(n,2) )<tol
            if abs( yf - node_coords(n,3) )<tol
                if abs( zf - node_coords(n,4) )<tol
                    Node_f = node_coords(n,1);
                end
            end
        end
        
    end
    
    NetList(j,:) = [j Node_i Node_m Zint(i)/2];
    NetList(j+1,:) = [j+1 0 Node_m Zg(i)];
    NetList(j+2,:) = [j+2 Node_m Node_f Zint(i)/2];
    
    j = j+3;
    
end

end