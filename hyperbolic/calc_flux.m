function [f_w,f_e] = calc_flux(fW,fP,fE,avrg,ge,gw)

  %Interpolate finite volume face values
  if strcmpi(avrg,'arithmetic')
    f_e = 0.5*(fE+fP);
    f_w = 0.5*(fW+fP);
  elseif strcmpi(avrg,'harmonic')
    f_e = (2.0*fE*fP)/(fE+fP);
    f_w = (2.0*fW*fP)/(fW+fP);
  elseif strcmpi(avrg,'geometric')
    f_e = ge*fE+(1-ge)*fP; %where: ge=(x_e-xP)/(xE-xP) 
    f_w = gw*fW+(1-gw)*fP; %where: gw=(x_w-xP)/(xW-xP)
  else
    error('Error: calc_flux (avrg)');
  end
