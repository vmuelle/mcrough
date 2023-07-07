function ippg(winkel, model)




model.MC.lightSource.sourceType   = 2; % 5: X/Y factorizable beam (e.g., a rectangular LED emitter)

model.MC.lightSource.theta     = 0; % [rad] Polar angle of direction the light collector is facing
model.MC.lightSource.phi       = 0; % [rad] Azimuthal angle of direction the light collector is facing

model.MC.useLightCollector        = true;
model.MC.lightCollector.x         = 0; % [cm] x position of either the center of the objective lens focal plane or the fiber tip
model.MC.lightCollector.y         = 0; % [cm] y position
model.MC.lightCollector.z         = 0; % [cm] z position

model.MC.lightCollector.theta     = (winkel*2*pi)/360; % [rad] Polar angle of direction the light collector is facing
model.MC.lightCollector.phi       = 0; % [rad] Azimuthal angle of direction the light collector is facing

model.MC.lightCollector.f         = 0.52; % [cm] Focal length of the objective lens (if light collector is a fiber, set this to Inf).
model.MC.lightCollector.diam      = 0.87; % [cm] Diameter of the light collector aperture. For an ideal thin lens, this is 2*f*tan(asin(NA)).
model.MC.lightCollector.fieldSize = model.G.Lx; % [cm] Field Size of the imaging system (diameter of area in object plane that gets imaged). Only used for finite f.
%model.MC.lightCollector.NA        = 0.22; % [-] Fiber NA. Only used for infinite f.

model.MC.lightCollector.res       = model.G.nx; % X and Y resolution of light collector in pixels, only used for finite f

model.MC.depositionCriteria.onlyCollected = false;



model = runMonteCarlo(model);
model = plot(model,'MC');
fluence = model.MC.normalizedFluenceRate(:,:,:);
transmittance = model.MC.normalizedIrradiance_zpos;
reflectance = model.MC.normalizedIrradiance_zneg;
image = model.MC.lightCollector.image;
name = strcat('ippg',string(winkel),'.mat');
save(name,'fluence','transmittance','reflectance','image');
%MCmatlab.plotAzFz(model1);

end


