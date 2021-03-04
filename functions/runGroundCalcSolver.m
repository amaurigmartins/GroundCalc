function [out] = runGroundCalcSolver(fname, cond_radius, rho_top, rho_bottom, h_top, i_src, t_fault, rho_cover, h_cover, plot_surf, plot_curr, plot_touch, plot_step)
addpath('../functions')

try
    ti=tic;
    dx=.5;
    dy=dx;
    SB=5;
    h = waitbar(0,'Parsing input data...','Name','Progress...');
    coords_src = importdata(fname);
    coords_src = coords_src.data;
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
    if MAX_LEN < 2*cond_radius
        MAX_LEN=4*cond_radius;
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
    cond_radius = cond_radius*ones([size(coords_src,1),1]);
    waitbar(1,h,'Parsing input data... Done!')
    
    waitbar(0,h,'Subdividing conductors...')
    coords_src=subdivcond([coords_src cond_radius],MAX_LEN);
    R=zeros(size(coords_src,1));
    l=coords_src(:,10);
    cond_radius=coords_src(:,11);
    waitbar(1,h,'Subdividing conductors... Done!')
    
    for o=1:size(coords_src,1)
        for j=1:size(coords_src,1)
            msg = sprintf('Computing matrix coefficient terms... Segment %i / %i.',o,j);
            waitbar(o/size(coords_src,1),h,msg);
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
    
    msg = sprintf('Inverting coefficient matrix %i x %i...',o,j);
    waitbar(0,h,msg);
    lambda=R\ones([size(coords_src,1),1]);
    msg = sprintf('Inverting coefficient matrix %i x %i... Done!',o,j);
    waitbar(1,h,msg);
    
    waitbar(0,h,'Computing leakage current densities...');
    V = i_src / sum(lambda.*l);
    delta=lambda*V;
    Ii=delta.*l;
    Rg=V/i_src;
    waitbar(1,h,'Computing leakage current densities... Done!')
    
    fig=1;
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
            clinep([coords_src(i,1) coords_src(i,4)+offx],[coords_src(i,2) coords_src(i,5)+offy],[coords_src(i,3) coords_src(i,6)+offz],[delta(i) delta(i)]);
            hold on;
        end
        colorbar;
        colormap(jet);
        grid on;
        hold off;
        title('Leakage current distribution')
        fig=fig+1;
    end
    
    waitbar(0,h,'Computing surface potentials...')
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
                msg = sprintf('Computing surface potentials... Point (%i , %i) / Source %i.',i,j,k);
                waitbar(i/size(X,1),h,msg);
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
    waitbar(1,h,'Computing surface potentia(ls... Done!')
    
    if plot_surf
        figure(fig)
        surf(X, Y, Us,'EdgeColor', 'None', 'facecolor', 'interp')
        xlim([min(min(X)) max(max(X))]);
        ylim([min(min(Y)) max(max(Y))]);
        view(2)
        axis on
        colorbar
        colormap(jet);
        title('Surface GPR distribution');
        fig=fig+1;
    end
    
    if plot_touch
        figure(fig)
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
        title('Touch voltages distribution');
        fig=fig+1;
    end
    
    if plot_step
        figure(fig)
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
        title('Step voltages distribution');
    end
    
    close(h);
    tf=toc(ti);
    
    fprintf('\n\nCOMPUTATION SUMMARY:\n--------------------\n');
    fprintf('    Top layer resistivity:                   rho_top = %g ohm.m\n',rho_top);
    fprintf('    Bottom layer resistivity:                rho_bottom = %g ohm.m\n',rho_bottom);
    fprintf('    Top layer thickness:                     h = %g m\n',h_top);
    fprintf('    Reflection coefficient:                  k = %g\n',((rho_top-rho_bottom)/(rho_top+rho_bottom)));
    fprintf('    Energization current:                    I = %g + j%g = %g A |_ %g°\n',real(i_src),imag(i_src),abs(i_src),rad2deg(angle(i_src)));
    fprintf('    Ground grid GPR:                         V = %g + j%g = %g V |_ %g°\n',real(V),imag(V),abs(V),rad2deg(angle(V)));
    fprintf('    Ground impedance:                        Rg = %g + j%g = %g ohms |_ %g°\n',real(Rg),imag(Rg),abs(Rg),rad2deg(angle(Rg)));
    fprintf('    Maximum surface GPR:                     Us,max = %g V\n',max(max(Us)))
    fprintf('    Fault time:                              t = %g s\n',t_fault);
    fprintf('    Cover material resistivity:              rho_cov = %g ohm.m\n',rho_cover);
    fprintf('    Cover layer thickness:                   h_cov = %g m\n',h_cover);
    fprintf('    Maximum allowable touch voltage:         Et,lim = %g V\n',Et_max);
    fprintf('    Maximum touch voltage:                   Et,max = %g V\n',max(max(Et)));
    fprintf('    Maximum allowable step voltage:          Ep,lim = %g V\n',Ep_max);
    fprintf('    Maximum step voltage:                    Ep,max = %g V\n',max(max(TotalGrad)));
    fprintf('    Number of source subdivisions:           Nf = %i \n',size(coords_src,1));
    fprintf('    Number of surface observation points:    Ns = %i \n\n', size(coordX_prof,2) * size(coordY_prof,2));
    disp(['    *** ELAPSED TIME:                      ' secs2hms(tf)]);
    
    save('groundcalc_data.mat');
    
catch MExc
    close(h);
end

end

