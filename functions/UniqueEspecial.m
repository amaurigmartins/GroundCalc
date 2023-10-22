function [ coords_src, cond_radius, cond_rho, cond_mur ] = UniqueEspecial( coords_src, cond_radius, cond_rho, cond_mur  )

%This function erases the identical segmentes from "coords_src" and the
%radius of the segment in "cond_radius".

for i=1:size(coords_src,1)
    if sum( abs(coords_src(i,:) - [0 0 0 0 0 0]) ) > 0.0001
        for j=(i+1):size(coords_src,1)

            if abs(coords_src(i,:) - coords_src(j,:)) < 0.0001
                coords_src(j,:) = 0;
                cond_radius(j) = 0;
                cond_rho(j) = 0;
                cond_mur(j) = 0;
            end

        end
    end
end

coords_src( all(~coords_src,2), : ) = [];
cond_radius( all(~cond_radius,2), : ) = [];
cond_rho( all(~cond_rho,2), : ) = [];
cond_mur( all(~cond_mur,2), : ) = [];

end

