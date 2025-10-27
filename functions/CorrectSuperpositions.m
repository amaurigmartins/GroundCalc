function [coords_src, cond_radius] = CorrectSuperpositions(coords_src, cond_radius)

% This function corrects every superposition found in the data inserted for
% conductors displacement.

warning('off','all');

if size(coords_src,1)>1                           % Only correct superpositions if there are multiple conductors

    % First erase the perfectly overlaid segments.
   
    [coords_src, cond_radius] = UniqueEspecial( coords_src, cond_radius );
    
    if size(coords_src,1)>1                       % Only correct superpositions if there still are multiple conductors

        % Deal with the other types of superpositions.

        [coords_src,cond_radius] = SolveOverlaid( coords_src, cond_radius );
    
    end

end

warning('on','all');

end