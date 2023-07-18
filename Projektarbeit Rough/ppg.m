function ppg(winkel, model,roughness_type)




model.MC.lightSource.sourceType   = 5; % 5: X/Y factorizable beam (e.g., a rectangular LED emitter)

model.MC.lightSource.focalPlaneIntensityDistribution.XDistr = 0; % X focal plane intensity distribution - 0: Top-hat
model.MC.lightSource.focalPlaneIntensityDistribution.XWidth = .02; % X focal plane in cm, 1/e^2 radius if top-hat

model.MC.lightSource.focalPlaneIntensityDistribution.YDistr = 0; % Y focal plane intensity distribution - 0: Top-hat
model.MC.lightSource.focalPlaneIntensityDistribution.YWidth = .01; % Y focal plane in cm, 1/e^2 radius if top-hat

model.MC.lightSource.angularIntensityDistribution.XDistr = 2; % X angular intensity distribution - 2: Cosine (Lambertian)
model.MC.lightSource.angularIntensityDistribution.YDistr = 2; % Y angular intensity distribution - 2: Cosine (Lambertian)


model.MC.lightSource.theta     = (winkel*2*pi)/360; % [rad] Polar angle of direction the light collector is facing
model.MC.lightSource.phi       = 0; % [rad] Azimuthal angle of direction the light collector is facing

model.MC.useLightCollector        = true;
model.MC.lightCollector.x         = 0; % [cm] x position of either the center of the objective lens focal plane or the fiber tip
model.MC.lightCollector.y         = 0.17; % [cm] y position
model.MC.lightCollector.z         = 0; % [cm] z position

model.MC.lightCollector.theta     = 0; % [rad] Polar angle of direction the light collector is facing
model.MC.lightCollector.phi       = 0; % [rad] Azimuthal angle of direction the light collector is facing

model.MC.lightCollector.f         = 1e-6; % [cm] Focal length of the objective lens (if light collector is a fiber, set this to Inf).
model.MC.lightCollector.diam      = 0.136; % [cm] Diameter of the light collector aperture. For an ideal thin lens, this is 2*f*tan(asin(NA)).
model.MC.lightCollector.fieldSize = 0.136;%model.G.Lx; % [cm] Field Size of the imaging system (diameter of area in object plane that gets imaged). Only used for finite f.
%model.MC.lightCollector.NA        = 0.22; % [-] Fiber NA. Only used for infinite f.

model.MC.lightCollector.res       = model.G.nx; % X and Y resolution of light collector in pixels, only used for finite f

model.MC.depositionCriteria.onlyCollected = true;



model = runMonteCarlo(model);
model = plot(model,'MC');
fluence = model.MC.normalizedFluenceRate(:,:,:);
transmittance = model.MC.normalizedIrradiance_zpos;
reflectance = model.MC.normalizedIrradiance_zneg;
image = model.MC.lightCollector.image;
name = strcat('outputs/ppg/',string(winkel),'deg',string(model.MC.wavelength),'um',string(roughness_type),'.mat');
save(name,'fluence','transmittance','reflectance','image');
%MCmatlab.plotAzFz(model1);

end


