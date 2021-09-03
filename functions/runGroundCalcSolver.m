function [out] = runGroundCalcSolver(coords_src, delta, buffer, cond_radius, rho_top, rho_bottom, h_top, i_src, t_fault, rho_cover, h_cover, plot_surf, plot_curr, plot_touch, plot_step, plot_grid, supress_messages, save_file)
% addpath('../functions')

try
    ti=tic;
    dx=delta;
    dy=dx;
    SB=buffer*delta;
    if ~ supress_messages
        h = waitbar(0,'Parsing input data...','Name','Progress...');
    end
%     coords_src = importdata(fname);
%     coords_src = coords_src.data;
    for i=1:size(coords_src,1)
        xi=coords_src(i,1);
        xf=coords_src(i,4);
        yi=coords_src(i,2);
        yf=coords_src(i,5);
        zi=coords_src(i,3);
        zf=coords_src(i,6);
        L(i) = sqrt((xf-xi)^2+(yf-yi)^2+(zf-zi)^2);
    end
    MAX_LEN = max(L) / 10;
    if MAX_LEN < 2*max(cond_radius)
        MAX_LEN=4*max(cond_radius);
    end
    minx=min(min(coords_src(:,1),coords_src(:,4)));
    maxx=max(max(coords_src(:,1),coords_src(:,4)));
    if minx==maxx
        coordX_prof=minx-SB:dx:maxx+SB;
    else
        coordX_prof=minx-SB:dx:maxx+SB;
    end
    miny=min(min(coords_src(:,2),coords_src(:,5)));
    maxy=max(max(coords_src(:,2),coords_src(:,5)));
    if miny==maxy
        coordY_prof=miny-SB:dy:maxy+SB;
    else
        coordY_prof=miny-SB:dy:maxy+SB;
    end
    t_fault=t_fault/1000;
%     cond_radius = cond_radius*ones([size(coords_src,1),1]);
    if ~ supress_messages
        waitbar(1,h,'Parsing input data... Done!')
    end
    
    if ~ supress_messages
        waitbar(0,h,'Subdividing conductors...')
    end
    coords_src=subdivcond([coords_src cond_radius],MAX_LEN);
    R=zeros(size(coords_src,1));
    l=coords_src(:,10);
    cond_radius=coords_src(:,11);
    if ~ supress_messages
        waitbar(1,h,'Subdividing conductors... Done!')
    end
    
    for o=1:size(coords_src,1)
        for j=1:size(coords_src,1)
            if ~ supress_messages
                msg = sprintf('Computing matrix coefficient terms... Segment %i / %i.',o,j);
                waitbar(o/size(coords_src,1),h,msg);
            end
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
                    this_z0=coords_src(j,9)+this_rad;
                else % vertical
                    this_x0=coords_src(j,7)+this_rad;
                    this_y0=coords_src(j,8);
                    this_z0=coords_src(j,9);
                end
            else
                this_x0=coords_src(j,7);
                this_y0=coords_src(j,8);
                this_z0=coords_src(j,9);
            end
            R(o,j)=greenfun2la(rho_top,rho_bottom,h_top,this_xs,this_ys,this_zs,this_xe,this_ye,this_ze,this_len,this_x0,this_y0,this_z0);
        end
    end
    
    if ~ supress_messages
        msg = sprintf('Inverting coefficient matrix %i x %i...',o,j);
        waitbar(0,h,msg);
    end
    lambda=R\ones([size(coords_src,1),1]);
    if ~ supress_messages
        msg = sprintf('Inverting coefficient matrix %i x %i... Done!',o,j);
        waitbar(1,h,msg);
    end
    if ~ supress_messages
        waitbar(0,h,'Computing leakage current densities...');
    end
    V = i_src / sum(lambda.*l);
    delta=lambda*V;
    Ii=delta.*l;
    Rg=V/i_src;
    if ~ supress_messages
        waitbar(1,h,'Computing leakage current densities... Done!')
    end
    if ~ supress_messages
        waitbar(0,h,'Computing surface potentials...')
    end
    [X,Y]=meshgrid(coordX_prof, coordY_prof);
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
                if ~ supress_messages
                    msg = sprintf('Computing surface potentials... Point (%i , %i) / Source %i.',i,j,k);
                    waitbar(i/size(X,1),h,msg);
                end
                this_x0=X(i,j);
                this_y0=Y(i,j);
                this_z0=0;
                this_xs=coords_src(k,1);
                this_ys=coords_src(k,2);
                this_zs=coords_src(k,3);
                this_xe=coords_src(k,4);
                this_ye=coords_src(k,5);
                this_ze=coords_src(k,6);
                %             this_rad=cond_radius(k,1);
                this_len=coords_src(k,10);
                GPR(i,j)= GPR(i,j)+delta(k)*greenfun2la(rho_top,rho_bottom,h_top,this_xs,this_ys,this_zs,this_xe,this_ye,this_ze,this_len,this_x0,this_y0,this_z0);
            end
        end
    end
    
    Us=abs(GPR);
    Et=abs(GPR-V);
    [Usx,Usy]=gradient(Us,1,1);
    TotalGrad=(Usx.^2+Usy.^2).^0.5;
    if ~ supress_messages
        waitbar(1,h,'Computing surface potentials... Done!')
    end
    
    if plot_surf
        figure(1)
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
        title('Surface GPR distribution');
    end
    
    if plot_touch
        figure(2)
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
        vmin=Et_max;
        vmax=max(max(Et));
        if vmax<vmin
            vmax=2*vmin;
        end
        caxis([vmin vmax])
        title(sprintf('Touch voltages distribution, E_{t,lim} = %1.0f V',Et_max));
    end
    
    if plot_step
        figure(3)
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
        vmin=Ep_max;
        vmax=max(max(TotalGrad));
        if vmax<vmin
            vmax=2*vmin;
        end
        caxis([vmin vmax])
        title(sprintf('Step voltages distribution, E_{p,lim} = %1.0f V',Ep_max));
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
            figure(4)
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
    end
    
    if  ~ supress_messages
        close(h);
    end
    
    tf=toc(ti);
    
    xx=[coords_src(:,1);coords_src(:,4)];
    yy=[coords_src(:,2);coords_src(:,5)];
    k = boundary(xx,yy);
    xv=xx(k);
    yv=yy(k);
    [in,on] = inpolygon(X,Y,xv,yv); % to find values inside the grid
    
     
    fprintf('\n\nCOMPUTATION SUMMARY:\n--------------------\n');
    fprintf('    Top layer resistivity:                   rho_top = %g ohm.m\n',rho_top);
    fprintf('    Bottom layer resistivity:                rho_bottom = %g ohm.m\n',rho_bottom);
    fprintf('    Top layer thickness:                     h = %g m\n',h_top);
    fprintf('    Reflection coefficient:                  k = %g\n',((rho_top-rho_bottom)/(rho_top+rho_bottom)));
    fprintf('    Energization current:                    I = %g + j%g = %g A %c %g°\n',real(i_src),imag(i_src),abs(i_src),char( hex2dec('2220') ),rad2deg(angle(i_src)));
    fprintf('    Ground grid GPR:                         V = %g + j%g = %g V %c %g°\n',real(V),imag(V),abs(V),char( hex2dec('2220') ),rad2deg(angle(V)));
    fprintf('    Ground impedance:                        Rg = %g + j%g = %g ohms %c %g°\n',real(Rg),imag(Rg),abs(Rg),char( hex2dec('2220') ),rad2deg(angle(Rg)));
    fprintf('    Maximum surface GPR:                     Us,max = %g V\n',max(max(Us)))
    fprintf('    Fault time:                              t = %g s\n',t_fault);
    fprintf('    Cover material resistivity:              rho_cov = %g ohm.m\n',rho_cover);
    fprintf('    Cover layer thickness:                   h_cov = %g m\n',h_cover);
    fprintf('    Maximum allowable touch voltage:         Et,lim = %g V\n',Et_max);
    fprintf('    Maximum touch voltage inside grid:       Et,max = %g V\n',max(Et(in|on)));
    fprintf('    Maximum allowable step voltage:          Ep,lim = %g V\n',Ep_max);
    fprintf('    Maximum step voltage inside grid:        Ep,max = %g V\n',max(TotalGrad(in|on)));
    fprintf('    Number of source subdivisions:           Nf = %i \n',size(coords_src,1));
    fprintf('    Number of surface observation points:    Ns = %i \n\n', size(coordX_prof,2) * size(coordY_prof,2));
    disp(['    *** ELAPSED TIME:                      ' secs2hms(tf)]);
    
    if save_file
        save('groundcalc_data.mat');
    end

    
    f = msgbox('All done! Check the Command Window for the computation summary.', 'GroundCalc','modal');
    
catch MExc
    close(h);
    end

end

