function [new_coords_src, new_cond_radius] = SolveOverlaid( coords_src, cond_radius )

% This function fix the imperfectly overlaid segments 
% in matrix coords_src. 

tol = 1E-5;

j = 1;

over = [0 0 0 0 0 0 0];

for n=1:size(coords_src,1)

    Pi_n = [coords_src(n,1) coords_src(n,2) coords_src(n,3)]; %Segment n initial point
    
    Pf_n = [coords_src(n,4) coords_src(n,5) coords_src(n,6)]; %Segment n final point
    
    for k=n+1:size(coords_src,1)  
            
        Pi_k = [coords_src(k,1) coords_src(k,2) coords_src(k,3)]; %Segment k initial point

        Pf_k = [coords_src(k,4) coords_src(k,5) coords_src(k,6)]; %Segment k final point

        check_Pin = PinSeg( Pi_n(1), Pi_n(2), Pi_n(3), Pi_k(1), Pi_k(2), Pi_k(3), Pf_k(1), Pf_k(2), Pf_k(3) ); %Check if initial point of segment n is in segment k

        check_Pfn = PinSeg( Pf_n(1), Pf_n(2), Pf_n(3), Pi_k(1), Pi_k(2), Pi_k(3), Pf_k(1), Pf_k(2), Pf_k(3) ); %Check if final point of segment n is in segment k
        
        check_Pik = PinSeg( Pi_k(1), Pi_k(2), Pi_k(3), Pi_n(1), Pi_n(2), Pi_n(3), Pf_n(1), Pf_n(2), Pf_n(3) ); %Check if initial point of segment n is in segment k

        check_Pfk = PinSeg( Pf_k(1), Pf_k(2), Pf_k(3), Pi_n(1), Pi_n(2), Pi_n(3), Pf_n(1), Pf_n(2), Pf_n(3) ); %Check if final point of segment n is in segment k

        
        if check_Pin                             % If initial point of segment n is in segment k

            if sum( abs(Pi_n - Pi_k) ) < tol    % Pi_n is at the beginning of segment k 
                over(j,:) = [0 n 0 k 0 0 0];
                j=j+1;
            elseif sum( abs(Pi_n - Pf_k) ) < tol % Pi_n is at the end of segment k 
                over(j,:) = [0 n 1 k 0 0 0];
                j=j+1;
            else % Pi_n is in segment k
                over(j,:) = [0 n 0.5 k 0 0 0];
                j=j+1;
            end

        end
        
        if check_Pik                             % If initial point of segment k is in segment n

            if sum( abs(Pi_k - Pi_n) ) < tol    % Pi_k is at the beginning of segment n 
                over(j,:) = [0 n 0 k 0 0 0];
                j=j+1;
            elseif sum( abs(Pi_k - Pf_n) ) < tol % Pi_k is at the end of segment n 
                over(j,:) = [1 n 0 k 0 0 0];
                j=j+1;
            else % Pi_k is in segment n
                over(j,:) = [0 k 0.5 n 0 0 0];
                j=j+1;
            end

        end

        if check_Pfn                          % If final point of segment n is in segment k

            if sum( abs(Pf_n - Pi_k) ) < tol % Pf_n is at the beginning of segment k 
                over(j,:) = [1 n 0 k 0 0 0];
                j=j+1;
            elseif sum( abs(Pf_n - Pf_k) ) < tol % Pf_n is at the end of segment k 
                over(j,:) = [1 n 1 k 0 0 0];
                j=j+1;
            else % Pf_n is in segment k
                over(j,:) = [1 n 0.5 k 0 0 0];
                j=j+1;
            end

        end
        
        if check_Pfk                          % If final point of segment k is in segment n

            if sum( abs(Pf_k - Pi_n) ) < tol % Pf_k is at the beginning of segment n 
                over(j,:) = [0 n 1 k 0 0 0];
                j=j+1;
            elseif sum( abs(Pf_k - Pf_n) ) < tol % Pf_k is at the end of segment n 
                over(j,:) = [1 n 1 k 0 0 0];
                j=j+1;
            else % Pf_k is in segment n
                over(j,:) = [1 k 0.5 n 0 0 0];
                j=j+1;
            end

        end
        
        if (~check_Pin) && (~check_Pfn) && (~check_Pik) && (~check_Pfk)  % If the initial or final point of n is not in k, check if they make a "X"
            
            Pis = [Pi_n(1) Pi_n(2) Pi_n(3); Pi_k(1) Pi_k(2) Pi_k(3)];
            Pfs = [Pf_n(1) Pf_n(2) Pf_n(3); Pf_k(1) Pf_k(2) Pf_k(3)];

            X = lineXline(Pis,Pfs); % Returns the intersection point of the lines, if they are parallel return NaN.

            X_in_n = PinSeg( X(1), X(2), X(3), Pi_n(1), Pi_n(2), Pi_n(3), Pf_n(1), Pf_n(2), Pf_n(3) ); % Checks if the intersection point is in the n segment
            X_in_k = PinSeg( X(1), X(2), X(3), Pi_k(1), Pi_k(2), Pi_k(3), Pf_k(1), Pf_k(2), Pf_k(3) ); % Checks if the intersection point is in the n segment
        
            
            if X_in_n && X_in_k % Segments n and k cross each other
                over(j,:) = [0.5 n 0.5 k X(1) X(2) X(3)];
                j=j+1;                
            end
        
        end
        
    end
    
end

over = unique(over, 'rows');

%Matrix "over"
%1st column: indicates if its the begginig (0), middle (0.5) or end (1) of segment n 
%2nd column: n
%3rd column: indicates if its the begginig (0), middle (0.5) or end (1) of segment k
%4th column: k
%5,6 and 7th column: point of intersection in case of "X" superpositions

new_coords_src = coords_src;     % Creates new matrix for the segments without overlaping
new_cond_radius = cond_radius;   % Creates new matrix for the segments without overlaping

if  sum( abs(over) ) < tol % No superpositions found
    
else

    for i=1:size(over,1)

    %     if n or k lines in coords_src = 0
    %         
    %         Faz nada
    %         
    %     else resolve aí

        if over(i,1) == 0 %Pi_n is in segment k

            % Where in segment k is the point Pi_n?

            if over(i,3) == 0   %Pi_n == Pi_k

                n = over(i,2);
                k = over(i,4);

                Pf_n = [coords_src(n,4) coords_src(n,5) coords_src(n,6)]; %Segment n final point

                Pi_k = [coords_src(k,1) coords_src(k,2) coords_src(k,3)]; %Segment k initial point

                Pf_k = [coords_src(k,4) coords_src(k,5) coords_src(k,6)]; %Segment k final point

                check_Pf = PinSeg( Pf_n(1), Pf_n(2), Pf_n(3), Pi_k(1), Pi_k(2), Pi_k(3), Pf_k(1), Pf_k(2), Pf_k(3) ); %Check if final point of segmente n is in segment k

                if check_Pf % If Pi_n and Pf_n are in segment k, we can erase segment n

                    new_coords_src(n,:) = 0;
                    new_cond_radius(n) = 0;

                end

            elseif over(i,3) == 0.5   %Pi_n in segment k 

                n = over(i,2);
                k = over(i,4);

                Pi_n = [coords_src(n,1) coords_src(n,2) coords_src(n,3)]; %Segment n initial point

                Pf_n = [coords_src(n,4) coords_src(n,5) coords_src(n,6)]; %Segment n final point

                Pi_k = [coords_src(k,1) coords_src(k,2) coords_src(k,3)]; %Segment k initial point

                Pf_k = [coords_src(k,4) coords_src(k,5) coords_src(k,6)]; %Segment k final point

                check_Pf = PinSeg( Pf_n(1), Pf_n(2), Pf_n(3), Pi_k(1), Pi_k(2), Pi_k(3), Pf_k(1), Pf_k(2), Pf_k(3) ); %Check if final point of segmente n is in segment k

                if check_Pf % If Pi_n and Pf_n are in segment k, we can erase segment n

                    new_coords_src(n,:) = 0;
                    new_cond_radius(n) = 0;

                else % Pi over Pi OR Pi over Pf OR "T" type intersection

                    check_Pi = PinSeg( Pi_k(1), Pi_k(2), Pi_k(3), Pi_n(1), Pi_n(2), Pi_n(3), Pf_n(1), Pf_n(2), Pf_n(3) ); %Check if final point of segmente n is in segment k
                    check_Pf = PinSeg( Pf_k(1), Pf_k(2), Pf_k(3), Pi_n(1), Pi_n(2), Pi_n(3), Pf_n(1), Pf_n(2), Pf_n(3) ); %Check if final point of segmente n is in segment k

                    if check_Pi % Pi over Pi

                        new_coords_src(n,:) = [Pi_k(1) Pi_k(2) Pi_k(3) Pf_n(1) Pf_n(2) Pf_n(3)];
                        new_coords_src(k,:) = [Pi_n(1) Pi_n(2) Pi_n(3) Pf_k(1) Pf_k(2) Pf_k(3)];

                        new_segment = [Pi_n(1) Pi_n(2) Pi_n(3) Pi_k(1) Pi_k(2) Pi_k(3)];
                        
                        if cond_radius(k)<cond_radius(n)    % Being conservative and choosing the smaller radius
                            new_radius = cond_radius(k);
                        else
                            new_radius = cond_radius(n);
                        end
                            
                        new_coords_src = [new_coords_src; new_segment];
                        new_cond_radius = [new_cond_radius; new_radius];

                    elseif check_Pf % Pi over Pf

                        new_coords_src(n,:) = [Pf_k(1) Pf_k(2) Pf_k(3) Pf_n(1) Pf_n(2) Pf_n(3)];
                        new_coords_src(k,:) = [Pi_k(1) Pi_k(2) Pi_k(3) Pi_n(1) Pi_n(2) Pi_n(3)];

                        new_segment = [Pi_n(1) Pi_n(2) Pi_n(3) Pf_k(1) Pf_k(2) Pf_k(3)];
                        
                        if cond_radius(k)<cond_radius(n)    % Being conservative and choosing the smaller radius
                            new_radius = cond_radius(k);
                        else
                            new_radius = cond_radius(n);
                        end 

                        new_coords_src = [new_coords_src; new_segment];
                        new_cond_radius = [new_cond_radius; new_radius];

                    else % "T" intersection

                        new_coords_src(k,:) = [Pi_k(1) Pi_k(2) Pi_k(3) Pi_n(1) Pi_n(2) Pi_n(3)];

                        new_segment = [Pi_n(1) Pi_n(2) Pi_n(3) Pf_k(1) Pf_k(2) Pf_k(3)];
                        new_radius = cond_radius(k); 

                        new_coords_src = [new_coords_src; new_segment];
                        new_cond_radius = [new_cond_radius; new_radius];

                    end

                end

            elseif over(i,3) == 1   %Pi_n == Pf_k

                n = over(i,2);
                k = over(i,4);

                Pf_n = [coords_src(n,4) coords_src(n,5) coords_src(n,6)]; %Segment n final point

                Pi_k = [coords_src(k,1) coords_src(k,2) coords_src(k,3)]; %Segment k initial point

                Pf_k = [coords_src(k,4) coords_src(k,5) coords_src(k,6)]; %Segment k final point

                check_Pf = PinSeg( Pf_n(1), Pf_n(2), Pf_n(3), Pi_k(1), Pi_k(2), Pi_k(3), Pf_k(1), Pf_k(2), Pf_k(3) ); %Check if final point of segmente n is in segment k

                if check_Pf % If Pi_n and Pf_n are in segment k, we can erase segment n

                    new_coords_src(n,:) = 0;
                    new_cond_radius(n) = 0;

                end

            end

        elseif over(i,1) == 1 %Pf_n is in segment k  

            % Where in segment k is the point Pf_n?

            if over(i,3) == 0   %Pf_n == Pi_k

                n = over(i,2);
                k = over(i,4);

                Pi_n = [coords_src(n,1) coords_src(n,2) coords_src(n,3)]; %Segment n initial point

                Pi_k = [coords_src(k,1) coords_src(k,2) coords_src(k,3)]; %Segment k initial point

                Pf_k = [coords_src(k,4) coords_src(k,5) coords_src(k,6)]; %Segment k final point

                check_Pi = PinSeg( Pi_n(1), Pi_n(2), Pi_n(3), Pi_k(1), Pi_k(2), Pi_k(3), Pf_k(1), Pf_k(2), Pf_k(3) ); %Check if final point of segmente n is in segment k

                if check_Pi % If Pi_n and Pf_n are in segment k, we can erase segment n

                    new_coords_src(n,:) = 0;
                    new_cond_radius(n) = 0;

                end

            elseif over(i,3) == 0.5   %Pf_n in segment k

                n = over(i,2);
                k = over(i,4);

                Pi_n = [coords_src(n,1) coords_src(n,2) coords_src(n,3)]; %Segment n initial point

                Pf_n = [coords_src(n,4) coords_src(n,5) coords_src(n,6)]; %Segment n final point

                Pi_k = [coords_src(k,1) coords_src(k,2) coords_src(k,3)]; %Segment k initial point

                Pf_k = [coords_src(k,4) coords_src(k,5) coords_src(k,6)]; %Segment k final point

                check_Pi = PinSeg( Pi_n(1), Pi_n(2), Pi_n(3), Pi_k(1), Pi_k(2), Pi_k(3), Pf_k(1), Pf_k(2), Pf_k(3) ); %Check if final point of segmente n is in segment k

                if check_Pi % If Pi_n and Pf_n are in segment k, we can erase segment n

                    new_coords_src(n,:) = 0;
                    new_cond_radius(n) = 0;

                else % Pf over Pi OR Pf over Pf OR "T" type intersection

                    check_Pi = PinSeg( Pi_k(1), Pi_k(2), Pi_k(3), Pi_n(1), Pi_n(2), Pi_n(3), Pf_n(1), Pf_n(2), Pf_n(3) ); %Check if final point of segmente n is in segment k
                    check_Pf = PinSeg( Pf_k(1), Pf_k(2), Pf_k(3), Pi_n(1), Pi_n(2), Pi_n(3), Pf_n(1), Pf_n(2), Pf_n(3) ); %Check if final point of segmente n is in segment k

                    if check_Pi % Pf over Pi

                        new_coords_src(n,:) = [Pi_n(1) Pi_n(2) Pi_n(3) Pi_k(1) Pi_k(2) Pi_k(3)];
                        new_coords_src(k,:) = [Pf_n(1) Pf_n(2) Pf_n(3) Pf_k(1) Pf_k(2) Pf_k(3)];

                        new_segment = [Pi_k(1) Pi_k(2) Pi_k(3) Pf_n(1) Pf_n(2) Pf_n(3)];
                        if cond_radius(k)<cond_radius(n)    % Being conservative and choosing the smaller radius
                            new_radius = cond_radius(k);
                        else
                            new_radius = cond_radius(n);
                        end 

                        new_coords_src = [new_coords_src; new_segment];
                        new_cond_radius = [new_cond_radius; new_radius];

                    elseif check_Pf % Pi over Pf

                        new_coords_src(n,:) = [Pf_k(1) Pf_k(2) Pf_k(3) Pi_n(1) Pi_n(2) Pi_n(3)];
                        new_coords_src(k,:) = [Pi_k(1) Pi_k(2) Pi_k(3) Pf_n(1) Pf_n(2) Pf_n(3)];

                        new_segment = [Pf_n(1) Pf_n(2) Pf_n(3) Pf_k(1) Pf_k(2) Pf_k(3)];
                        if cond_radius(k)<cond_radius(n)    % Being conservative and choosing the smaller radius
                            new_radius = cond_radius(k);
                        else
                            new_radius = cond_radius(n);
                        end

                        new_coords_src = [new_coords_src; new_segment];
                        new_cond_radius = [new_cond_radius; new_radius];

                    else % "T" intersection

                        new_coords_src(k,:) = [Pi_k(1) Pi_k(2) Pi_k(3) Pf_n(1) Pf_n(2) Pf_n(3)];

                        new_segment = [Pf_n(1) Pf_n(2) Pf_n(3) Pf_k(1) Pf_k(2) Pf_k(3)];
                        new_radius = cond_radius(k); 

                        new_coords_src = [new_coords_src; new_segment];
                        new_cond_radius = [new_cond_radius; new_radius];

                    end

                end

            elseif over(i,3) == 1   %Pf_n == Pf_k

                n = over(i,2);
                k = over(i,4);

                Pi_n = [coords_src(n,1) coords_src(n,2) coords_src(n,3)]; %Segment n final point

                Pi_k = [coords_src(k,1) coords_src(k,2) coords_src(k,3)]; %Segment k initial point

                Pf_k = [coords_src(k,4) coords_src(k,5) coords_src(k,6)]; %Segment k final point

                check_Pi = PinSeg( Pi_n(1), Pi_n(2), Pi_n(3), Pi_k(1), Pi_k(2), Pi_k(3), Pf_k(1), Pf_k(2), Pf_k(3) ); %Check if final point of segmente n is in segment k

                if check_Pi % If Pi_n and Pf_n are in segment k, we can erase segment n

                    new_coords_src(n,:) = 0;
                    new_cond_radius(n) = 0;

                end

            end

        elseif over(i,1) == 0.5 %Segments n and k make a "X"

            n = over(i,2);
            k = over(i,4);

            Pi_n = [coords_src(n,1) coords_src(n,2) coords_src(n,3)]; %Segment n initial point
            Pf_n = [coords_src(n,4) coords_src(n,5) coords_src(n,6)]; %Segment n final point
            Pi_k = [coords_src(k,1) coords_src(k,2) coords_src(k,3)]; %Segment n initial point
            Pf_k = [coords_src(k,4) coords_src(k,5) coords_src(k,6)]; %Segment n final point

            X = [over(i,5) over(i,6) over(i,7)]; %Intersection point

            new_coords_src(n,:) = [Pi_n(1) Pi_n(2) Pi_n(3) X(1) X(2) X(3)];
            new_coords_src(k,:) = [Pi_k(1) Pi_k(2) Pi_k(3) X(1) X(2) X(3)];

            new_segment = [X(1) X(2) X(3) Pf_n(1) Pf_n(2) Pf_n(3);
                           X(1) X(2) X(3) Pf_k(1) Pf_k(2) Pf_k(3)];
            new_radius = [cond_radius(n);
                          cond_radius(k)];

            new_coords_src = [new_coords_src;new_segment];
            new_cond_radius = [new_cond_radius;new_radius];

        end

    end

end
    
new_coords_src( all(~new_coords_src,2), : ) = [];
new_cond_radius( all(~new_cond_radius,2), : ) = [];

[new_coords_src, new_cond_radius] = UniqueEspecial( new_coords_src, new_cond_radius );

end