function [ out ] = greenfunwrapper( rho1, rho2, h, xs, ys, zs, xe, ye, ze, comp, x0, y0, z0, enforce_segmentationonly )

 out = greenfun2la(rho1,rho2,h,xs,ys,zs,xe,ye,ze,comp,x0,y0,z0,enforce_segmentationonly);

end

