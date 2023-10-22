function [ out ] = greenfunwrapper( rho1, rho2, h, xs, ys, zs, xe, ye, ze, comp, x0, y0, z0, enforce_segmentationonly )

if enforce_segmentationonly
    out=[];
else
    out = greenfun2la_(rho1,rho2,h,xs,ys,zs,xe,ye,ze,comp,x0,y0,z0);
end

end

