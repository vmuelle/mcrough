close all
clear all
clc

%% Geometry definition
MCmatlab.closeMCmatlabFigures();
model = MCmatlab.model;

model.G.nx                = 200; % Number of bins in the x direction
model.G.ny                = 200; % Number of bins in the y direction
model.G.nz                = 200; % Number of bins in the z direction
model.G.Lx                = 1; % x size of simulation cuboid in cm
model.G.Ly                = 1; % y size of simulation cuboid in cm
model.G.Lz                = 0.822; % z size of simulation cuboid in cm

model.G.mediaPropertiesFunc = @mediaPropertiesFunc; % Media properties defined as a function at the end of this file
model.G.geomFunc          = @geometryDefinition; % Function to use for defining the distribution of media in the cuboid.

model = plot(model,'G');

%% Monte Carlo simulation
model.MC.simulationTimeRequested  = .1; % Time duration of the simulation in min
model.MC.nExamplePaths            = 100;
model.MC.matchedInterfaces        = false;
model.MC.boundaryType             = 1; % 1: All cuboid boundaries are escaping
model.MC.wavelength               = 940; % Excitation wavelength in nm, used for determination of optical properties for excitation light

model.MC.lightSource.sourceType   = 5; % 5: X/Y factorizable beam (e.g., a rectangular LED emitter)

model.MC.lightSource.focalPlaneIntensityDistribution.XDistr = 0; % X focal plane intensity distribution - 0: Top-hat
model.MC.lightSource.focalPlaneIntensityDistribution.XWidth = .02; % X focal plane in cm, 1/e^2 radius if top-hat

model.MC.lightSource.focalPlaneIntensityDistribution.YDistr = 0; % Y focal plane intensity distribution - 0: Top-hat
model.MC.lightSource.focalPlaneIntensityDistribution.YWidth = .01; % Y focal plane in cm, 1/e^2 radius if top-hat

model.MC.lightSource.angularIntensityDistribution.XDistr = 2; % X angular intensity distribution - 2: Cosine (Lambertian)
model.MC.lightSource.angularIntensityDistribution.YDistr = 2; % Y angular intensity distribution - 2: Cosine (Lambertian)

model = runMonteCarlo(model);
model = plot(model,'MC');
IR = model.MC.normalizedFluenceRate(:,:,:);
save('IR.mat','IR');

%% Geometry function(s) (see readme for details)
% A geometry function takes as input X,Y,Z matrices as returned by the
% "ndgrid" MATLAB function as well as any parameters the user may have
% provided in the definition of Ginput. It returns the media matrix M,
% containing numerical values indicating the media type (as defined in
% mediaPropertiesFunc) at each voxel location.
function M = geometryDefinition(X,Y,Z,parameters)
  zsurf = 0.005;
  SC_thick = 0.002;
  LE_thick = 0.025;
  PD_thick = 0.015;
  RD_thick = 0.15;
  Hypo_thick = 0.1;
  Muskel_thick = 0.225;
  M = ones(size(X)); % air
  M(Z > zsurf) = 2; % stratum corneum
  M(Z > zsurf + SC_thick) = 3; % living epidermis
  M(Z > zsurf + SC_thick + LE_thick) = 4; % stratum papillare
  M(Z > zsurf + SC_thick + LE_thick + PD_thick) = 5; % stratum reticulare
  M(Z > zsurf + SC_thick + LE_thick + PD_thick + RD_thick) = 6; % hypodermis
  M(Z > zsurf + SC_thick + LE_thick + PD_thick + RD_thick + Hypo_thick) = 7; % muskel
  M(Z > zsurf + SC_thick + LE_thick + PD_thick + RD_thick + Hypo_thick + Muskel_thick) = 8; % knochen
end

%% Media Properties function (see readme for details)
% The media properties function defines all the optical and thermal
% properties of the media involved by filling out and returning a
% "mediaProperties" array of "mediumProperties" objects with various
% properties. The j indices are the numbers that are referred to in the
% geometry function (in this case, 1 for "air" and 2 for "standard tissue")
% See the readme file or the examples for a list of properties you may
% specify. Most properties may be specified as a numeric constant or as
% function handles.
% 
% The function must take one input; the cell array containing any special
% parameters you might specify above in the model file, for example
% parameters that you might loop over in a for loop. In most simulations
% this "parameters" cell array is empty. Dependence on wavelength is shown
% in examples 4 and 23. Dependence on excitation fluence rate FR,
% temperature T or fractional heat damage FD can be specified as in
% examples 12-15.
function mediaProperties = mediaPropertiesFunc(parameters)
  mediaProperties = MCmatlab.mediumProperties;

  j=1;
  mediaProperties(j).name  = 'Air';
  mediaProperties(j).mua   = 1e-8; % Absorption coefficient in cm^-1
  mediaProperties(j).mus   = 1e-8; % Scattering coefficient in cm^-1
  mediaProperties(j).g     = 1; % Henyey-Greenstein scattering anisotropy
  mediaProperties(j).n     = 1; % Refractive index

  j=2;
  lambda = 940; % wavelength in nm
  SC_Vw = 0.05; % fraction of water
  mua_w_940 = 0.36; % absorption coefficient for water at 940 nm in cm^-1
  mua_baseline = 7.84*10^7*lambda^-3.255; 
  g_SC = 0.942;
  mediaProperties(j).name  = 'Stratum Corneum';
% mediaProperties(j).mua   = ((0.1-0.3*10^-4*lambda)+0.125*mua_baseline)*(1-SC_Vw)*mua_w_940; % Absorption coefficient in cm^-1 
  mediaProperties(j).mua   = 0.8; % Absorption coefficient in cm^-1
  mediaProperties(j).mus   = 33.6/(1-g_SC); % Scattering coefficient in cm^-1
  mediaProperties(j).g     = 0.942; % Henyey-Greenstein scattering anisotropy
  mediaProperties(j).n     = 1.5; % Refractive index

  j=3;
  LE_Vmel = 0.11; % fraction of melanin
  LE_Vw = 0.2; % fraction of water
  mua_mel_940 = 83; % absorption coefficient for melanin at 940 nm in cm^-1
  g_LE = 0.8;
  mediaProperties(j).name  = 'Living Epidermis';
% mediaProperties(j).mua   = (LE_Vmel*mua_mel_940+(1-LE_Vmel)*mua_baseline)*(1-LE_Vw)+LE_Vw*mua_w_940; % Absorption coefficient in cm^-1
  mediaProperties(j).mua   = 0.8; % Absorption coefficient in cm^-1
  mediaProperties(j).mus   = 33.6/(1-g_LE); % Scattering coefficient in cm^-1
  mediaProperties(j).g     = 0.8; % Henyey-Greenstein scattering anisotropy
  mediaProperties(j).n     = 1.34; % Refractive index

  j=4;
  g_PD = 0.9;
  mediaProperties(j).name  = 'Stratum Papillare';
  mediaProperties(j).mua   = 0.51; % Absorption coefficient in cm^-1
  mediaProperties(j).mus   = 15.1/(1-g_PD); % Scattering coefficient in cm^-1
  mediaProperties(j).g     = 0.9; % Henyey-Greenstein scattering anisotropy
  mediaProperties(j).n     = 1.4; % Refractive index

  j=5;
  g_RD = 0.8;
  mediaProperties(j).name  = 'Stratum Reticulare';
  mediaProperties(j).mua   = 0.55; % Absorption coefficient in cm^-1
  mediaProperties(j).mus   = 15.1/(1-g_RD); % Scattering coefficient in cm^-1
  mediaProperties(j).g     = 0.8; % Henyey-Greenstein scattering anisotropy
  mediaProperties(j).n     = 1.4; % Refractive index

  j=6;
  mediaProperties(j).name  = 'Hypodermis';
  mediaProperties(j).mua   = 0.91; % Absorption coefficient in cm^-1
  mediaProperties(j).mus   = 109.9; % Scattering coefficient in cm^-1
  mediaProperties(j).g     = 0.75; % Henyey-Greenstein scattering anisotropy
  mediaProperties(j).n     = 1.44; % Refractive index

  j=7;
  mediaProperties(j).name  = 'Muskel';
  mediaProperties(j).mua   = 1.61; % Absorption coefficient in cm^-1
  mediaProperties(j).mus   = 79.8; % Scattering coefficient in cm^-1
  mediaProperties(j).g     = 0.9112; % Henyey-Greenstein scattering anisotropy
  mediaProperties(j).n     = 1.41; % Refractive index

  j=8;
  mediaProperties(j).name  = 'Knochen';
  mediaProperties(j).mua   = 0.244; % Absorption coefficient in cm^-1
  mediaProperties(j).mus   = 215.7; % Scattering coefficient in cm^-1
  mediaProperties(j).g     = 0.93; % Henyey-Greenstein scattering anisotropy
  mediaProperties(j).n     = 1.41; % Refractive index
end