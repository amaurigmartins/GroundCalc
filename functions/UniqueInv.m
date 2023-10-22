function [coords_src, cond_radius, cond_rho, cond_mur ] = UniqueInv(coords_src, cond_radius, cond_rho, cond_mur )

%This function erases the identical inverted segmentes from "coords_src" and the
%radius of the segment in "cond_radius".

for n=1:size(coords_src,1)

    if sum( abs(coords_src(n,:) - [0 0 0 0 0 0]) ) > 0.0001
    
        Pi_n = [coords_src(n,1) coords_src(n,2) coords_src(n,3)]; %Segment n initial point

        Pf_n = [coords_src(n,4) coords_src(n,5) coords_src(n,6)]; %Segment n final point

        for k=(n+1):size(coords_src,1)  

                Pi_k = [coords_src(k,1) coords_src(k,2) coords_src(k,3)]; %Segment k initial point

                Pf_k = [coords_src(k,4) coords_src(k,5) coords_src(k,6)]; %Segment k final point

                if (abs(Pi_n - Pf_k) < 0.0001)
                    if (abs(Pf_n - Pi_k) < 0.0001)
                        coords_src(k,:) = 0;
                        cond_radius(k) = 0;
                        cond_rho(k) = 0;
                        cond_mur(k) = 0;
                    end
                end

        end
    
    end
    
end

coords_src( all(~coords_src,2), : ) = [];
cond_radius( all(~cond_radius,2), : ) = [];
cond_rho( all(~cond_rho,2), : ) = [];
cond_mur( all(~cond_mur,2), : ) = [];

end