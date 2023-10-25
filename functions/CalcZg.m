function [Zg] = CalcZg(coords_src, rho_top, rho_bottom, h_top, supress_messages, enforce_segmentationonly, app)

% This function uses the matrix method to calculate the impedance to
% the ground from each segment in the grounding grid.
timestamp = @(t) datetime(t,'ConvertFrom','datenum'); %because reasons

R=zeros(size(coords_src,1)); %Cria a matriz R e a preenche com zeros
l=coords_src(:,10); %Cria um vetor com os comprimentos de cada segmento
cond_radius=coords_src(:,11); %Cria um vetor com o raio de cada segmento

if ~ supress_messages
    app.ProgrammessagesTextArea.Value{end+1} = [sprintf('%s: Computing matrix coefficient terms. Hang on...',timestamp(now))];
end

for o=1:size(coords_src,1)      %Loop que preenche a matriz R
    for j=1:size(coords_src,1)

        this_xs=coords_src(o,1);
        this_ys=coords_src(o,2);
        this_zs=coords_src(o,3);
        this_xe=coords_src(o,4);
        this_ye=coords_src(o,5);
        this_ze=coords_src(o,6);
        this_rad=cond_radius(o,1);
        this_len=coords_src(o,10);

        if o==j
            if abs(this_zs-this_ze) <= 0.0001 % horizontal
                this_x0=coords_src(j,7);
                this_y0=coords_src(j,8);
                this_z0=coords_src(j,9)+this_rad;%Caso o segmento seja horizontal o z0 é posicionado na superfície de cima dele
            else % vertical
                this_x0=coords_src(j,7)+this_rad;%Caso o segmento seja vertical o x0 é posicionado na superfície dele
                this_y0=coords_src(j,8);
                this_z0=coords_src(j,9);
            end
        else
            this_x0=coords_src(j,7); %x0 y0 z0:ponto médio do segmento
            this_y0=coords_src(j,8);
            this_z0=coords_src(j,9);
        end
        
        R(o,j)=greenfunwrapper(rho_top,rho_bottom,h_top,this_xs,this_ys,this_zs,this_xe,this_ye,this_ze,this_len,this_x0,this_y0,this_z0,enforce_segmentationonly);
%       O elemento R_o,j é a própria função de Greens.
    end
    if ~ supress_messages
        app.ProgrammessagesTextArea.Value{end+1} = [sprintf('-- Terms %d to %d... Done!',((o-1)*j)+1,o*j)];
    end
end

%comments are for the weak
R=enforce_symm(R);

if ~ supress_messages
    app.ProgrammessagesTextArea.Value{end+1} = [sprintf('%s: Inverting coefficients matrix...',timestamp(now))];
end
%Eq. (6.22)
lambda=R\ones([size(coords_src,1),1]); %Multiplica o inverso da matriz R por um vetor de 1s. É a mesma coisa que inv(R)*ones([size(coords_src,1),1])

Zg = zeros( size(lambda,1), 1 );

if ~ supress_messages
    app.ProgrammessagesTextArea.Value{end+1} = [sprintf('%s: Computing lossless leakage current densities and earth admittances...',timestamp(now))];
end

for i=1:size(lambda,1)
   Yg = lambda(i)*l(i); %Fazendo isso vou ter Ij/GPR, que é a admitância pro solo.
   Zg(i) = 1/Yg;
end

end