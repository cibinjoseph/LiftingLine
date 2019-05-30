function liftNonDimTheoretical = getLiftnonDimTheoretical(AR,cl_slope,semispan,y)
% returns theoretical value of non-dim lift of an elliptical wing
% sectional lift/(0.5*rho*Vinf^2*mean chord*cl_slope*theta)

liftNonDimTheoretical = 4*AR/(pi*AR+cl_slope)*sqrt(1-(y/semispan).^2);
return
