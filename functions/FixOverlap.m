function [new_coords_src, new_cond_radius, new_cond_rho, new_cond_mur, NumOver] = FixOverlap( coords_src, cond_radius, cond_rho, cond_mur )

orig_state = warning;warning('off','all');


NumOver = 0;

%This function fix the imperfectly and perfectly overlaid segments
%in matrix coords_src.

tol = 0.0001;

%First whe erase the perfectly overlaid segments.

[coords_src, cond_radius, cond_rho, cond_mur ] = UniqueEspecial( coords_src, cond_radius, cond_rho, cond_mur  );

[coords_src, cond_radius, cond_rho, cond_mur ] = UniqueInv(coords_src, cond_radius, cond_rho, cond_mur );

%Then we deal with the other types of superpositions.

j=1;

for n=1:size(coords_src,1)

    Pi_n = [coords_src(n,1) coords_src(n,2) coords_src(n,3)]; %Segment n initial point

    Pf_n = [coords_src(n,4) coords_src(n,5) coords_src(n,6)]; %Segment n final point

    for k=1:size(coords_src,1)
        if k ~= n

            Pi_k = [coords_src(k,1) coords_src(k,2) coords_src(k,3)]; %Segment k initial point

            Pf_k = [coords_src(k,4) coords_src(k,5) coords_src(k,6)]; %Segment k final point


            check_Pi = PinSeg( Pi_n(1), Pi_n(2), Pi_n(3), Pi_k(1), Pi_k(2), Pi_k(3), Pf_k(1), Pf_k(2), Pf_k(3) ); %Check if initial point of segmente n is in segment k

            check_Pf = PinSeg( Pf_n(1), Pf_n(2), Pf_n(3), Pi_k(1), Pi_k(2), Pi_k(3), Pf_k(1), Pf_k(2), Pf_k(3) ); %Check if final point of segmente n is in segment k

            Pis = [Pi_n(1) Pi_n(2) Pi_n(3); Pi_k(1) Pi_k(2) Pi_k(3)];
            Pfs = [Pf_n(1) Pf_n(2) Pf_n(3); Pf_k(1) Pf_k(2) Pf_k(3)];

            X = lineXline(Pis,Pfs);

            X_in_n = PinSeg( X(1), X(2), X(3), Pi_n(1), Pi_n(2), Pi_n(3), Pf_n(1), Pf_n(2), Pf_n(3) );
            X_in_k = PinSeg( X(1), X(2), X(3), Pi_k(1), Pi_k(2), Pi_k(3), Pf_k(1), Pf_k(2), Pf_k(3) );

            if check_Pi

                if abs(Pi_n - Pi_k) < tol %Pi_n is at the beginning of segment k
                    %                     fprintf('O ponto inicial do segmento %i está no ponto inicial do segmento %i.\n', n, k);
                    over(j,:) = [0 n 0 k 0 0 0];
                    j=j+1;
                elseif abs(Pi_n - Pf_k) < tol %Pi_n is at the end of segment k
                    %                     fprintf('O ponto inicial do segmento %i está no ponto final do segmento %i.\n', n, k);
                    over(j,:) = [0 n 1 k 0 0 0];
                    j=j+1;
                else %Pi_n is in segment k
                    %                     fprintf('O ponto inicial do segmento %i está no meio do segmento %i.\n', n, k);
                    over(j,:) = [0 n 0.5 k 0 0 0];
                    j=j+1;
                end

            elseif check_Pf

                if abs(Pf_n - Pi_k) < tol %Pf_n is at the beginning of segment k
                    %                     fprintf('O ponto final do segmento %i está no ponto inicial do segmento %i.\n', n, k);
                    over(j,:) = [1 n 0 k 0 0 0];
                    j=j+1;
                elseif abs(Pf_n - Pf_k) < tol %Pf_n is at the end of segment k
                    %                     fprintf('O ponto final do segmento %i está no ponto final do segmento %i.\n', n, k);
                    over(j,:) = [1 n 1 k 0 0 0];
                    j=j+1;
                else %Pf_n is in segment k
                    %                     fprintf('O ponto final do segmento %i está no meio do segmento %i.\n', n, k);
                    over(j,:) = [1 n 0.5 k 0 0 0];
                    j=j+1;
                end

            elseif X_in_n && X_in_k

                %                 fprintf('Os segmentos %i e %i se cruzam no ponto %f %f %f.\n', n, k, X(1), X(2), X(3));
                over(j,:) = [0.5 n 0.5 k X(1) X(2) X(3)];
                j=j+1;

            end

        end
    end

end

%Matrix "over"
%1st column: indicates if its the begginig (0), middle (0.5) or end (1) of segment n
%2nd column: n
%3rd column: indicates if its the begginig (0), middle (0.5) or end (1) of segment k
%4th column: k
%5,6 and 7th column: point of intersection in case of "X" superpositions

% coords_src

% over

new_coords_src = coords_src;
new_cond_radius = cond_radius;
new_cond_rho = cond_rho;
new_cond_mur = cond_mur;

for i=1:size(over,1)

    if over(i,1) == 0 %Pi_n is in segment k
        if over(i,3) == 0   %Pi_n == Pi_k

            n = over(i,2);
            k = over(i,4);

            Pf_n = [coords_src(n,4) coords_src(n,5) coords_src(n,6)]; %Segment n final point

            Pi_k = [coords_src(k,1) coords_src(k,2) coords_src(k,3)]; %Segment k initial point

            Pf_k = [coords_src(k,4) coords_src(k,5) coords_src(k,6)]; %Segment k final point

            check_Pf = PinSeg( Pf_n(1), Pf_n(2), Pf_n(3), Pi_k(1), Pi_k(2), Pi_k(3), Pf_k(1), Pf_k(2), Pf_k(3) ); %Check if final point of segmente n is in segment k

            if check_Pf
                new_coords_src(n,:) = 0;
                new_cond_radius(n) = 0;
                new_cond_rho(n) = 0;
                new_cond_mur(n) = 0;
            end

        end

        if over(i,3) == 0.5   %Pi_n in segment k

            NumOver = NumOver+1;

            n = over(i,2);
            k = over(i,4);

            Pi_n = [coords_src(n,1) coords_src(n,2) coords_src(n,3)]; %Segment n initial point

            Pf_n = [coords_src(n,4) coords_src(n,5) coords_src(n,6)]; %Segment n final point

            Pi_k = [coords_src(k,1) coords_src(k,2) coords_src(k,3)]; %Segment k initial point

            Pf_k = [coords_src(k,4) coords_src(k,5) coords_src(k,6)]; %Segment k final point

            check_Pf = PinSeg( Pf_n(1), Pf_n(2), Pf_n(3), Pi_k(1), Pi_k(2), Pi_k(3), Pf_k(1), Pf_k(2), Pf_k(3) ); %Check if final point of segmente n is in segment k

            if check_Pf
                new_coords_src(n,:) = 0;
                new_cond_radius(n) = 0;
                new_cond_rho(n) = 0;
                new_cond_mur(n) = 0;
            else
                new_coords_src(k,:) = [Pi_k(1) Pi_k(2) Pi_k(3) Pi_n(1) Pi_n(2) Pi_n(3)];

                new_segment = [Pi_n(1) Pi_n(2) Pi_n(3) Pf_k(1) Pf_k(2) Pf_k(3)];

                Check_Pi_K = PinSeg( Pi_k(1), Pi_k(2), Pi_k(3), Pi_n(1), Pi_n(2), Pi_n(3), Pf_n(1), Pf_n(2), Pf_n(3) );
                Check_Pf_K = PinSeg( Pf_k(1), Pf_k(2), Pf_k(3), Pi_n(1), Pi_n(2), Pi_n(3), Pf_n(1), Pf_n(2), Pf_n(3) );

                if Check_Pi_K
                    if cond_radius(k) < cond_radius(n) %It's a "Pi over Pi" type intersection
                        new_cond_radius(k) = cond_radius(k);
                        new_cond_rho(k) = cond_rho(k);
                        new_cond_mur(k) = cond_mur(k);
                        new_radius = cond_radius(k);
                        new_rho = cond_rho(k);
                        new_mur = cond_mur(k);
                    else
                        new_cond_radius(k) = cond_radius(n);
                        new_cond_rho(k) = cond_rho(n);
                        new_cond_mur(k) = cond_mur(n);
                        new_radius = cond_radius(k);
                        new_rho = cond_rho(k);
                        new_mur = cond_mur(k);
                    end
                elseif Check_Pf_K
                    if cond_radius(k) < cond_radius(n) %It's a "Pi over Pf" type intersection
                        new_cond_radius(k) = cond_radius(k);
                        new_cond_rho(k) = cond_rho(k);
                        new_cond_mur(k) = cond_mur(k);
                        new_radius = cond_radius(k);
                        new_rho = cond_rho(k);
                        new_mur = cond_mur(k);
                    else
                        new_cond_radius(k) = cond_radius(k);
                        new_cond_rho(k) = cond_rho(k);
                        new_cond_mur(k) = cond_mur(k);
                        new_radius = cond_radius(n);
                        new_rho = cond_rho(n);
                        new_mur = cond_mur(n);
                    end
                else
                    new_radius = cond_radius(k); %It's a "T" type intersection
                    new_rho = cond_rho(k); %It's a "T" type intersection
                    new_mur = cond_mur(k); %It's a "T" type intersection
                end

                new_coords_src = [new_coords_src;new_segment];
                new_cond_radius = [new_cond_radius;new_radius];
                new_cond_rho = [new_cond_rho;new_rho];
                new_cond_mur = [new_cond_mur;new_mur];
            end

        end

        if over(i,3) == 1   %Pi_n == Pf_k

            n = over(i,2);
            k = over(i,4);

            Pf_n = [coords_src(n,4) coords_src(n,5) coords_src(n,6)]; %Segment n final point

            Pi_k = [coords_src(k,1) coords_src(k,2) coords_src(k,3)]; %Segment k initial point

            Pf_k = [coords_src(k,4) coords_src(k,5) coords_src(k,6)]; %Segment k final point

            check_Pf = PinSeg( Pf_n(1), Pf_n(2), Pf_n(3), Pi_k(1), Pi_k(2), Pi_k(3), Pf_k(1), Pf_k(2), Pf_k(3) ); %Check if final point of segmente n is in segment k

            if check_Pf
                new_coords_src(n,:) = 0;
                new_cond_radius(n) = 0;
                new_cond_rho(n) = 0;
                new_cond_mur(n) = 0;
            end

        end
    elseif over(i,1) == 1 %Pf_n is in segment k
        if over(i,3) == 0   %Pf_n == Pi_k

            n = over(i,2);
            k = over(i,4);

            Pi_n = [coords_src(n,1) coords_src(n,2) coords_src(n,3)]; %Segment n initial point

            Pi_k = [coords_src(k,1) coords_src(k,2) coords_src(k,3)]; %Segment k initial point

            Pf_k = [coords_src(k,4) coords_src(k,5) coords_src(k,6)]; %Segment k final point

            check_Pi = PinSeg( Pi_n(1), Pi_n(2), Pi_n(3), Pi_k(1), Pi_k(2), Pi_k(3), Pf_k(1), Pf_k(2), Pf_k(3) ); %Check if final point of segmente n is in segment k

            if check_Pi
                new_coords_src(n,:) = 0;
                new_cond_radius(n) = 0;
                new_cond_rho(n) = 0;
                new_cond_mur(n) = 0;
            end

        end

        if over(i,3) == 0.5   %Pf_n in segment k

            NumOver = NumOver+1;

            n = over(i,2);
            k = over(i,4);

            Pi_n = [coords_src(n,1) coords_src(n,2) coords_src(n,3)]; %Segment n initial point

            Pf_n = [coords_src(n,4) coords_src(n,5) coords_src(n,6)]; %Segment n final point

            Pi_k = [coords_src(k,1) coords_src(k,2) coords_src(k,3)]; %Segment k initial point

            Pf_k = [coords_src(k,4) coords_src(k,5) coords_src(k,6)]; %Segment k final point

            check_Pi = PinSeg( Pi_n(1), Pi_n(2), Pi_n(3), Pi_k(1), Pi_k(2), Pi_k(3), Pf_k(1), Pf_k(2), Pf_k(3) ); %Check if final point of segmente n is in segment k

            if check_Pi
                new_coords_src(n,:) = 0;
                new_cond_radius(n) = 0;
                new_cond_rho(n) = 0;
                new_cond_mur(n) = 0;
            else
                new_coords_src(k,:) = [Pi_k(1) Pi_k(2) Pi_k(3) Pf_n(1) Pf_n(2) Pf_n(3)];

                new_segment = [Pf_n(1) Pf_n(2) Pf_n(3) Pf_k(1) Pf_k(2) Pf_k(3)];

                Check_Pi_K = PinSeg( Pi_k(1), Pi_k(2), Pi_k(3), Pi_n(1), Pi_n(2), Pi_n(3), Pf_n(1), Pf_n(2), Pf_n(3) );
                Check_Pf_K = PinSeg( Pf_k(1), Pf_k(2), Pf_k(3), Pi_n(1), Pi_n(2), Pi_n(3), Pf_n(1), Pf_n(2), Pf_n(3) );

                if Check_Pi_K
                    if cond_radius(k) < cond_radius(n) %It's a "Pf over Pi" type intersection
                        new_cond_radius(k) = cond_radius(k);
                        new_cond_rho(k) = cond_rho(k);
                        new_cond_mur(k) = cond_mur(k);
                        new_radius = cond_radius(k);
                        new_rho = cond_rho(k);
                        new_mur = cond_mur(k);
                    else
                        new_cond_radius(k) = cond_radius(n);
                        new_cond_rho(k) = cond_rho(n);
                        new_cond_mur(k) = cond_mur(n);
                        new_radius = cond_radius(k);
                        new_rho = cond_rho(k);
                        new_mur = cond_mur(k);
                    end
                elseif Check_Pf_K
                    if cond_radius(k) < cond_radius(n) %It's a "Pf over Pf" type intersection
                        new_cond_radius(k) = cond_radius(k);
                        new_cond_rho(k) = cond_rho(k);
                        new_cond_mur(k) = cond_mur(k);
                        new_radius = cond_radius(k);
                        new_rho = cond_rho(k);
                        new_mur = cond_mur(k);
                    else
                        new_cond_radius(k) = cond_radius(k);
                        new_cond_rho(k) = cond_rho(k);
                        new_cond_mur(k) = cond_mur(k);
                        new_radius = cond_radius(n);
                        new_rho = cond_rho(n);
                        new_mur = cond_mur(n);
                    end
                else
                    new_radius = cond_radius(k); %It's a "T" type intersection
                    new_rho = cond_rho(k); %It's a "T" type intersection
                    new_mur = cond_mur(k); %It's a "T" type intersection
                end

                new_coords_src = [new_coords_src;new_segment];
                new_cond_radius = [new_cond_radius;new_radius];
                new_cond_rho = [new_cond_rho;new_rho];
                new_cond_mur = [new_cond_mur;new_mur];
            end

        end

        if over(i,3) == 1   %Pf_n == Pf_k

            n = over(i,2);
            k = over(i,4);

            Pi_n = [coords_src(n,1) coords_src(n,2) coords_src(n,3)]; %Segment n final point

            Pi_k = [coords_src(k,1) coords_src(k,2) coords_src(k,3)]; %Segment k initial point

            Pf_k = [coords_src(k,4) coords_src(k,5) coords_src(k,6)]; %Segment k final point

            check_Pi = PinSeg( Pi_n(1), Pi_n(2), Pi_n(3), Pi_k(1), Pi_k(2), Pi_k(3), Pf_k(1), Pf_k(2), Pf_k(3) ); %Check if final point of segmente n is in segment k

            if check_Pi
                new_coords_src(n,:) = 0;
                new_cond_radius(n) = 0;
                new_cond_rho(n) = 0;
                new_cond_mur(n) = 0;
            end

        end
    elseif over(i,1) == 0.5 %Segments n and k make a "X"

        NumOver = NumOver+1;

        n = over(i,2);

        Pi_n = [coords_src(n,1) coords_src(n,2) coords_src(n,3)]; %Segment n initial point

        Pf_n = [coords_src(n,4) coords_src(n,5) coords_src(n,6)]; %Segment n final point

        X = [over(i,5) over(i,6) over(i,7)]; %Intersection point

        new_coords_src(n,:) = [Pi_n(1) Pi_n(2) Pi_n(3) X(1) X(2) X(3)];

        new_segment = [X(1) X(2) X(3) Pf_n(1) Pf_n(2) Pf_n(3)];
        new_radius = cond_radius(n);
        new_rho = cond_rho(n);
        new_mur = cond_mur(n);

        new_coords_src = [new_coords_src;new_segment];
        new_cond_radius = [new_cond_radius;new_radius];
        new_cond_rho = [new_cond_rho;new_rho];
        new_cond_mur = [new_cond_mur;new_mur];

    end

end

% NumOver

new_coords_src( all(~new_coords_src,2), : ) = [];
new_cond_radius( all(~new_cond_radius,2), : ) = [];
new_cond_rho( all(~new_cond_rho,2), : ) = [];
new_cond_mur( all(~new_cond_mur,2), : ) = [];

[new_coords_src, new_cond_radius, new_cond_rho, new_cond_mur] = UniqueEspecial( new_coords_src, new_cond_radius, new_cond_rho, new_cond_mur );

[new_coords_src, new_cond_radius, new_cond_rho, new_cond_mur] = UniqueInv( new_coords_src, new_cond_radius, new_cond_rho, new_cond_mur );

warning(orig_state);


end