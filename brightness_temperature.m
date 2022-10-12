%% This function computes brightness temperturature


% By Andrew John Buggee
%%

function Tb = brightness_temperature(radiance, spectral_variable, spectral_units)

con = physical_constants();

if strcmp(spectral_units, 'microns')==true

    spec_var = spectral_variable./1e6;                      % meters - converted from microns
    L = radiance.*1e6;                                      % W/m^2/sr/meter - converted to a spectral dependence using meters

    Tb = con.h * con.c./(spec_var * con.k_B) * 1/(log(2*con.h * con.c^2/(spec_var.^5 .*L) +1));

else

    error([newline,'I dont recognize the spectral units of your input.',newline])

end



end