function [out] = runGroundCalcSolver(fname, cond_radius, rho_top, rho_bottom, h_top, i_src, t_fault, rho_cover, h_cover, plot_surf, plot_curr, plot_touch, plot_step)

ti=tic;

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
    coordX_prof=minx-10:1:maxx+10;
else
    coordX_prof=minx-abs(maxx-minx)*.1:1:maxx+abs(maxx-minx)*.1;
end
miny=min(min(coords_src(:,2),coords_src(:,5)));
maxy=max(max(coords_src(:,2),coords_src(:,5)));
if miny==maxy
    coordY_prof=miny-10:1:maxy+10;
else
    coordY_prof=miny-abs(maxy-miny)*.1:1:maxy+abs(maxy-miny)*.1;
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
        this_rad=raiocond(o,1);
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
            R(o,j)=somagama(rho_top,rho_bottom,h_top,this_xs,this_ys,this_zs,this_xe,this_ye,this_ze,this_len,this_x0,this_y0,this_z0);
    end
end

msg = sprintf('Inverting coefficient matrix %i x %i...',o,j);
waitbar(0,h,msg);
lambda=R\ones([size(coords_src,1),1])
msg = sprintf('Inverting coefficient matrix %i x %i... Done!',o,j);
waitbar(1,h,msg);

waitbar(0,h,'Computing leakage current densities...');
V = i_src / sum(lambda.*l);
delta=lambda*V;
Ii=delta.*l;
Rg=V/I;
waitbar(1,h,'Computing leakage current densities... Done!')

fig=1;

if plot_curr
    figure(fig);
    for i=1:size(coords_src,1); color_line([coords_src(i,1) coords_src(i,4) ],[coords_src(i,2) coords_src(i,5) ],abs(delta(i)),'LineWidth',3); hold on; end;
    colorbar;
    colormap(jet);
    grid on;
    hold off;
    title('Leakage current distribution')
    fig=fig+1;
end
end

