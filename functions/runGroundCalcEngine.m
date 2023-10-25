if ~ supress_messages  
    app.ProgrammessagesTextArea.Value{end+1} = [sprintf('%s: Parsing input data...',timestamp(now))];
end

orig_state = warning;warning('off','all');
this_dir=fullfile(WORKDIR,jobid);
mkdir(this_dir);
warning(orig_state);
cd(this_dir);

TOL=1e-6;
t_fault=t_fault/1000; 
FaultNodeCoords =  coords_src(src_idx,1:3);

X = [coords_src(:, 1); coords_src(:, 4)]';
Y = [coords_src(:, 2); coords_src(:, 5)]';
Z = [coords_src(:, 3); coords_src(:, 6)]';

% Calculate bounds based on user inputs
if min(X) == max(X)  % Single vertical conductor
    x0 = mean(mean(X)) - buffer;
    x1 = mean(mean(X)) + buffer;
else
    x0 = min(min(X)) - buffer;
    x1 = max(max(X)) + buffer;
end

if min(Y) == max(Y)  % Single vertical conductor
    y0 = mean(mean(Y)) - buffer;
    y1 = mean(mean(Y)) + buffer;
else
    y0 = min(min(Y)) - buffer;
    y1 = max(max(Y)) + buffer;
end

% Generate the profile lines
dx=delta;
dy=dx;
coordX_prof = linspace(x0, x1, ceil((x1 - x0) / dx));
coordY_prof = linspace(y0, y1, ceil((y1 - y0) / dy));

% Gather X, Y, Z into a single 3 x N matrix
all_coords = [X; Y; Z];
% Replicate the FaultNodeCoords to make it the same size as all_coords
fault_matrix = repmat(FaultNodeCoords', 1, size(all_coords, 2));
% Calculate the distances in one fell swoop
distances = sqrt(sum((all_coords - fault_matrix).^2, 1));
% Find the index of the most distant point. Oh, how they grow up so fast.
[~, most_distant_idx] = max(distances);
% Extract the coordinates of the prodigal point
ExtremityNodeCoords = all_coords(:, most_distant_idx)';

if abs(i_ret)<=TOL
    ReturnNodeCoords = [0 0 0];              % Coordenadas da corrente de retorno pelo neutro
else
    ReturnNodeCoords = coords_src(ret_idx,1:3);              % Coordenadas da corrente de retorno pelo neutro
end

ti=tic;

previous_file=fullfile(this_dir,'grid_data.mat');
if isfile(previous_file)
    testdata=load(previous_file);


    if ~(((testdata.enforce_equipotential==enforce_equipotential) && (testdata.enforce_segmentationonly==enforce_segmentationonly)) && (isfield(testdata,'Zbus') && ~enforce_equipotential))
        enforce_rebuild=true;
        if ~ supress_messages
            app.ProgrammessagesTextArea.Value{end+1} = [sprintf('%s: Found a ''grid_data.mat'' file from a previous input session, but current computation settings are different. Let us rebuild this big boy...',timestamp(now))];
        end
    else
        if ~ supress_messages
            app.ProgrammessagesTextArea.Value{end+1} = [sprintf('%s: Found a ''grid_data.mat'' file from a previous input session. Let us take a shortcut...',timestamp(now))];
        end
    end
end

if isfile(previous_file) && ~enforce_rebuild
    load(previous_file);
else
    %    Fixing overlapping conductors
    if ~ supress_messages
        app.ProgrammessagesTextArea.Value{end+1} = [sprintf('%s: Checking for overlapping conductors. This may take a while...',timestamp(now))];
    end

    if size(coords_src,1)>1
        NumOver = 1; %Assumimos que existe pelo menos uma superposição
        i=0;

        while(NumOver) %A função é chamada até que nenhum segmento esteja sobre o outro
            [coords_src,cond_radius,cond_rho,cond_mur,NumOver] = FixOverlap( coords_src, cond_radius, cond_rho, cond_mur );
            i=i+1;
        end
    end
    coords_src = AlinhaCond(coords_src);

    if ~ supress_messages
        app.ProgrammessagesTextArea.Value{end+1} = [sprintf('%s: Subdividing conductors...',timestamp(now))];
    end

    L = zeros( size(coords_src,1), 1 );
    for i=1:size(coords_src,1)  % Calcula o tamanho de cada condutor
        xi=coords_src(i,1);
        xf=coords_src(i,4);
        yi=coords_src(i,2);
        yf=coords_src(i,5);
        zi=coords_src(i,3);
        zf=coords_src(i,6);
        L(i) = sqrt((xf-xi)^2+(yf-yi)^2+(zf-zi)^2);
    end

    if MAX_LEN < 2*max(cond_radius) % Essa  condição existe pois precisamos que o diametro do segmento 2*max(cond_radius) seja pequeno
        % se comparado ao seu comprimento.
        MAX_LEN=4*max(cond_radius);     % Com essa linha o tamanho máximo de cada segmento se torna o dobro do maior diametro da malha.
    end

    coords_src = subdivcond([coords_src cond_radius cond_rho cond_mur],MAX_LEN); %Realiza a subdivisão dos condutores
    %A função acima retorna uma matriz onde cada linha contem:
    %[xi yi zi xf yf zf xm ym zm L rad rho mur]
    %xi yi zi: coordenadas iniciais do segmento
    %xf yf zf: coordenadas finais do segmento
    %xm ym zm: coordenadas do ponto médio do segmento
    %L: comprimento do segmento
    %rad: raio do segmento

    %redefine radii, rho and mur for the subdivided grid
    cond_radius=coords_src(:,11);
    cond_rho=coords_src(:,12);
    cond_mur=coords_src(:,13);

    save(previous_file,'coords_src','cond_radius','cond_rho','cond_mur','enforce_equipotential','enforce_segmentationonly');

    if ~ supress_messages
        app.ProgrammessagesTextArea.Value{end+1} = [sprintf('%s: *** Progress saved into file ''grid_data.mat''. Future runs will recover from this point. You may redefine your system by ticking the box ''Rebuild circuit model from scratch'' in the ''Computations'' tab.',timestamp(now))];
    end
end

if ~enforce_equipotential

    if ~(isfile(previous_file) && ~enforce_rebuild)
        % Construindo a NetList
        node_coords = NodeCoordsCalc( coords_src );

        if ~ supress_messages
            app.ProgrammessagesTextArea.Value{end+1} = [sprintf('%s: Computing leakage impedances...',timestamp(now))];
        end
        Zg = CalcZg(coords_src, rho_top, rho_bottom, h_top, supress_messages, enforce_segmentationonly, app);

        if ~ supress_messages
            app.ProgrammessagesTextArea.Value{end+1} = [sprintf('%s: Computing internal impedances...',timestamp(now))];
        end
        Zint = CalcZint(coords_src,freq,cond_rho,cond_mur,cond_radius);

        if ~ supress_messages
            app.ProgrammessagesTextArea.Value{end+1} = [sprintf('%s: Building netlist, Ybus and Zbus matrices...',timestamp(now))];
        end
        NetList = CalcNetList(coords_src, node_coords, Zg, Zint);
        % Calculando Zbus
        NumBus = size(node_coords,1);
        Ybus = YbusFromNet(NetList, 0, NumBus);
        Zbus = inv(Ybus);

        if ~ supress_messages
            app.ProgrammessagesTextArea.Value{end+1} = [sprintf('%s: *** Progress saved into file ''grid_data.mat''. Future runs will recover from this point. You may redefine your system by ticking the box ''Rebuild circuit model from scratch'' in the ''Computations'' tab.',timestamp(now))];
        end
        save(previous_file,'node_coords','Zg','Zint','NetList','NumBus','Ybus','Zbus','-append');
    end

    % Calculando a tensão em cada nó
    if ~ supress_messages
        app.ProgrammessagesTextArea.Value{end+1} = [sprintf('%s: Computing node voltages...',timestamp(now))];
    end
    I = zeros(size(Zbus,1),1);

    FaultNode = FindNodeIndex(FaultNodeCoords(1),FaultNodeCoords(2),FaultNodeCoords(3),node_coords);

    I(FaultNode) = i_src;   % Aplicando a corrente no ponto desejado

    if abs(i_ret) > TOL
        ReturnNode = FindNodeIndex(ReturnNodeCoords(1),ReturnNodeCoords(2),ReturnNodeCoords(3),node_coords);
        I(ReturnNode) = -i_ret;   % Aplicando a corrente de retorno pelo neutro
    end

    V = Ybus\I;

    % Calculando as correntes dispersadas em cada segmento
    if ~ supress_messages
        app.ProgrammessagesTextArea.Value{end+1} = [sprintf('%s: Computing segment leakage currents...',timestamp(now))];
    end
    Ii = CalcLeakageIi(coords_src, node_coords, V, Zg);
    clear delta; delta = zeros( size(Ii,1) );

    for i=1:size(Ii,1)
        delta(i) = abs(Ii(i))/coords_src(i,10);
    end

    % Cálculo da resistência da malha de aterramento
    if ~ supress_messages
        app.ProgrammessagesTextArea.Value{end+1} = [sprintf('%s: Computing grounding impedance...',timestamp(now))];
    end
    Rg=V(FaultNode)/i_src; %Eq. (6.26). Cálculo da resistência da malha de aterramento

else
    if ~ supress_messages
        app.ProgrammessagesTextArea.Value{end+1} = [sprintf('%s: Running equipotential engine (BUT WHY?)...',timestamp(now))];
    end
    [V, delta, Ii, Rg] = CalcEqui(coords_src, i_src, rho_top, rho_bottom, h_top, supress_messages, enforce_segmentationonly, app);
    
end

% Gráficos e Resultados

    if ~ supress_messages
        app.ProgrammessagesTextArea.Value{end+1} = [sprintf('%s: Computing surface potentials...',timestamp(now))];
    end

[X,Y]=meshgrid(coordX_prof, coordY_prof);

% X é uma matriz com as coordenadas de minx-5 até maxx+5 espaçadas de 0.5
% em todas as linhas
% Y é uma matriz em que cada linha possui o mesmo número repetidas vezes
% sua primeira linha é preenchida por minx-5, a segunda minx-5+0.5, etc

Ib=0.116/sqrt(t_fault);

if(h_cover==0)
    Cs=1;
    R2Fs=6*Cs*rho_top;
    R2Fp=1.5*Cs*rho_top;
else
    Cs=1-(0.106*((1-(rho_top/rho_cover))/(2*h_cover+0.106)));
    R2Fs=6*Cs*rho_cover;
    R2Fp=1.5*Cs*rho_cover;
end

Rb=1000;
Ep_max=(Rb+R2Fs)*Ib;
Et_max=(Rb+R2Fp)*Ib;
GPR=zeros(size(X,1),size(X,2));

for i=1:size(X,1)
    for j=1:size(Y,2)
        for k=1:size(coords_src,1)
            this_x0=X(i,j);
            this_y0=Y(i,j);
            this_z0=0;
            this_xs=coords_src(k,1);
            this_ys=coords_src(k,2);
            this_zs=coords_src(k,3);
            this_xe=coords_src(k,4);
            this_ye=coords_src(k,5);
            this_ze=coords_src(k,6);
            this_len=coords_src(k,10);
            %         GPR(i,j)= GPR(i,j)+delta(k)*NewGreenFunction(rho_top,rho_bottom,h_top,this_xs,this_ys,this_zs,this_xe,this_ye,this_ze,this_len,this_x0,this_y0,this_z0);
            GPR(i,j)= GPR(i,j)+delta(k)*greenfunwrapper(rho_top,rho_bottom,h_top,this_xs,this_ys,this_zs,this_xe,this_ye,this_ze,this_len,this_x0,this_y0,this_z0,enforce_segmentationonly);
        end

    end
    if ~ supress_messages
        app.ProgrammessagesTextArea.Value{end+1} = [sprintf('-- Terms %d to %d... Done!',((i-1)*j*k)+1,i*j*k)];
    end
end

% Calculate the distance from each point in the grid to FaultNodeCoords
distances = sqrt((X - FaultNodeCoords(1)).^2 + (Y - FaultNodeCoords(2)).^2);
% Find the indices of the minimum distance
[~, min_idx] = min(distances(:));
% Translate the linear index back to 2D indices
[row_idx, col_idx] = ind2sub(size(distances), min_idx);
% Fetch the corresponding GPR value
FaultPointGPR = GPR(row_idx, col_idx);

% Calculate the distance from each point in the grid to ExtremityNodeCoords
distances = sqrt((X - ExtremityNodeCoords(1)).^2 + (Y - ExtremityNodeCoords(2)).^2);
% Find the indices of the minimum distance
[~, min_idx] = min(distances(:));
% Translate the linear index back to 2D indices
[row_idx, col_idx] = ind2sub(size(distances), min_idx);
% Fetch the corresponding GPR value
ExtremityPointGPR = GPR(row_idx, col_idx);

Us=abs(GPR);
Et=abs(GPR-min(V));  %Considering the worst case cenario
[Usx,Usy]=gradient(Us,1,1);
TotalGrad=(Usx.^2+Usy.^2).^0.5;

if plot_surf
    f=figure(1);
    if plot_grid
        for i=1:size(coords_src,1); plot3([coords_src(i,1) coords_src(i,4)],[coords_src(i,2) coords_src(i,5)],[max(max(Us))*1.05 max(max(Us))*1.05],'k-','LineWidth',1.5); hold on; end;
    end
    surf(X, Y, Us,'EdgeColor', 'None', 'facecolor', 'interp')
    xlim([min(min(X)) max(max(X))]);
    ylim([min(min(Y)) max(max(Y))]);
    view(2)
    axis on
    colorbar
    colormap(jet);
    title(sprintf('Surface GPR Distribution, GPR_{max} = %1.0f V',max(max(Us))));
    cameratoolbar(f,'show');
end

if plot_touch
    f=figure(2);
    if plot_grid
        for i=1:size(coords_src,1); plot3([coords_src(i,1) coords_src(i,4)],[coords_src(i,2) coords_src(i,5)],[max(max(Et))*1.05 max(max(Et))*1.05],'k-','LineWidth',1.5); hold on; end;
    end
    surf(X, Y, Et,'EdgeColor', 'None', 'facecolor', 'interp')
    xlim([min(min(X)) max(max(X))]);
    ylim([min(min(Y)) max(max(Y))]);
    view(2)
    axis on
    colorbar
    colormap(jet);
    vmin=min(min(Et));
    vmax=max(max(Et));
    if vmax<vmin
        vmax=2*vmin;
    end
    caxis([vmin vmax])
    title(sprintf('Touch voltages distribution, E_{t,lim} = %1.0f V',Et_max));
    cameratoolbar(f,'show');
end

if plot_step
    f=figure(3);
    if plot_grid
        for i=1:size(coords_src,1); plot3([coords_src(i,1) coords_src(i,4)],[coords_src(i,2) coords_src(i,5)],[max(max(TotalGrad))*1.05 max(max(TotalGrad))*1.05],'k-','LineWidth',1.5); hold on; end;
    end
    surf(X, Y, TotalGrad,'EdgeColor', 'None', 'facecolor', 'interp')
    xlim([min(min(X)) max(max(X))]);
    ylim([min(min(Y)) max(max(Y))]);
    view(2)
    axis on
    colorbar
    colormap(jet);
    vmin=min(min(TotalGrad));
    vmax=max(max(TotalGrad));
    if vmax<vmin
        vmax=2*vmin;
    end
    caxis([vmin vmax])
    title(sprintf('Step voltages distribution, E_{p,lim} = %1.0f V',Ep_max));
    cameratoolbar(f,'show');
end



if plot_curr
    for i=1:size(coords_src,1)
        if abs(coords_src(i,3)-coords_src(i,6))<0.0001
            offx=0;
            offy=0;
            offz=0;
        else
            offx=0.0001;
            offy=0.0002;
            offz=0.0003;
        end
        f=figure(4);
        clinep([coords_src(i,1) coords_src(i,4)+offx],[coords_src(i,2) coords_src(i,5)+offy],[coords_src(i,3) coords_src(i,6)+offz],[delta(i) delta(i)]);
        hold on;
    end
    view([45 45]);
    colorbar;
    colormap(jet);
    grid on;
    hold off;
    title('Leakage currents distribution')
    axis equal
    cameratoolbar(f,'show');
end

tf=toc(ti);

xx=[coords_src(:,1);coords_src(:,4)];
yy=[coords_src(:,2);coords_src(:,5)];
k = boundary(xx,yy);
xv=xx(k);
yv=yy(k);
[in,on] = inpolygon(X,Y,xv,yv); % to find values inside the grid

if enforce_equipotential;equi_str='YES';else; equi_str='NO';end
if enforce_segmentationonly;ptsrc_str='YES';else; ptsrc_str='NO';end

app.ProgrammessagesTextArea.Value{end+1} = [sprintf('%s: Process finished!',timestamp(now))];
app.ProgrammessagesTextArea.Value{end+1} = [sprintf('\n\nCOMPUTATION SUMMARY:\n----------------------------------------\n')];
app.ProgrammessagesTextArea.Value{end+1} = [sprintf('Equipotential assumption: %s \n', equi_str)];
app.ProgrammessagesTextArea.Value{end+1} = [sprintf('Point source model: %s \n', ptsrc_str)];
app.ProgrammessagesTextArea.Value{end+1} = [sprintf('Top layer resistivity: %g ohm.m\n', rho_top)];
app.ProgrammessagesTextArea.Value{end+1} = [sprintf('Bottom layer resistivity: %g ohm.m\n', rho_bottom)];
app.ProgrammessagesTextArea.Value{end+1} = [sprintf('Top layer thickness: %g m\n', h_top)];
app.ProgrammessagesTextArea.Value{end+1} = [sprintf('Reflection coefficient: %g\n', ((rho_top-rho_bottom)/(rho_top+rho_bottom)))];
app.ProgrammessagesTextArea.Value{end+1} = [sprintf('Energization current: %g + j%g = %g A %c %g°\n', real(i_src), imag(i_src), abs(i_src), char(hex2dec('2220')), rad2deg(angle(i_src)))];
app.ProgrammessagesTextArea.Value{end+1} = [sprintf('Maximum surface GPR: %g V\n', max(max(Us)))];
app.ProgrammessagesTextArea.Value{end+1} = [sprintf('Fault clearance time: %g s\n', t_fault)];
app.ProgrammessagesTextArea.Value{end+1} = [sprintf('Cover material resistivity: %g ohm.m\n', rho_cover)];
app.ProgrammessagesTextArea.Value{end+1} = [sprintf('Cover layer thickness: %g m\n', h_cover)];
app.ProgrammessagesTextArea.Value{end+1} = [sprintf('Maximum allowable touch voltage: %g V\n', Et_max)];
app.ProgrammessagesTextArea.Value{end+1} = [sprintf('Maximum touch voltage inside grid: %g V\n', max(Et(in|on)))];
app.ProgrammessagesTextArea.Value{end+1} = [sprintf('Maximum allowable step voltage: %g V\n', Ep_max)];
app.ProgrammessagesTextArea.Value{end+1} = [sprintf('Maximum step voltage inside grid: %g V\n', max(TotalGrad(in|on)))];
app.ProgrammessagesTextArea.Value{end+1} = [sprintf('Number of source subdivisions: %i \n', size(coords_src, 1))];
app.ProgrammessagesTextArea.Value{end+1} = [sprintf('Number of surface observation points: %i \n', size(coordX_prof, 2) * size(coordY_prof, 2))];


if ~enforce_equipotential

    fault = FindNodeIndex(FaultNodeCoords(1),FaultNodeCoords(2),FaultNodeCoords(3),node_coords);
    ext = FindNodeIndex(ExtremityNodeCoords(1), ExtremityNodeCoords(2), ExtremityNodeCoords(3), node_coords);

    app.ProgrammessagesTextArea.Value{end+1} = [sprintf('Grounding impedance: %g + j%g = %g ohms %c %g°\n', real(Rg), imag(Rg), abs(Rg), char(hex2dec('2220')), rad2deg(angle(Rg)))];
    app.ProgrammessagesTextArea.Value{end+1} = [sprintf('Voltage at fault point: %f %c %f\n', abs(V(fault)), char(hex2dec('2220')), rad2deg(angle(V(fault))))];
    app.ProgrammessagesTextArea.Value{end+1} = [sprintf('Voltage at grid extremity: %f %c %f\n', abs(V(ext)), char(hex2dec('2220')), rad2deg(angle(V(ext))))];
    app.ProgrammessagesTextArea.Value{end+1} = [sprintf('GPR at fault point: %f\n', FaultPointGPR)];
    app.ProgrammessagesTextArea.Value{end+1} = [sprintf('GPR at grid extremity:  %f\n', ExtremityPointGPR)];


else

    app.ProgrammessagesTextArea.Value{end+1} = [sprintf('Grounding impedance: %g + j%g = %g ohms %c %g°\n', real(Rg), imag(Rg), abs(Rg), char(hex2dec('2220')), rad2deg(angle(Rg)))];
    app.ProgrammessagesTextArea.Value{end+1} = [sprintf('Grounding grid GPR: %g + j%g = %g V %c %g°\n', real(V), imag(V), abs(V), char(hex2dec('2220')), rad2deg(angle(V)))];

end

app.ProgrammessagesTextArea.Value{end+1} = [sprintf('*** ELAPSED TIME: %s\n\n\n', secs2hms(tf))];

if save_mat
    save(fullfile(this_dir,'groundcalc_results.mat'));
    FolderName = fullfile(this_dir,'plots');
    mkdir(FolderName)
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    for iFig = 1:length(FigList)
        FigHandle = FigList(iFig);
        FigName   = get(FigHandle, 'Name');
        savefig(FigHandle, fullfile(FolderName,[FigName '.fig']));
    end
end


f = msgbox('All done! Check the ''Program messages'' window for the computation summary.', 'GroundCalc','modal');