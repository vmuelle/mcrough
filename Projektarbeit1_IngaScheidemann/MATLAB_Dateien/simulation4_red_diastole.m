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
model.MC.wavelength               = 700; % Excitation wavelength in nm, used for determination of optical properties for excitation light

model.MC.lightSource.sourceType   = 5; % 5: X/Y factorizable beam (e.g., a rectangular LED emitter)

model.MC.lightSource.focalPlaneIntensityDistribution.XDistr = 0; % X focal plane intensity distribution - 0: Top-hat
model.MC.lightSource.focalPlaneIntensityDistribution.XWidth = .02; % X focal plane in cm, 1/e^2 radius if top-hat

model.MC.lightSource.focalPlaneIntensityDistribution.YDistr = 0; % Y focal plane intensity distribution - 0: Top-hat
model.MC.lightSource.focalPlaneIntensityDistribution.YWidth = .01; % Y focal plane in cm, 1/e^2 radius if top-hat

model.MC.lightSource.angularIntensityDistribution.XDistr = 2; % X angular intensity distribution - 2: Cosine (Lambertian)
model.MC.lightSource.angularIntensityDistribution.YDistr = 2; % Y angular intensity distribution - 2: Cosine (Lambertian)


model = runMonteCarlo(model);
model = plot(model,'MC');
red_diastole = model.MC.normalizedFluenceRate(:,:,:);
save('red_diastole.mat','red_diastole');

%% Geometry function(s)
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

%% Media Properties function
function mediaProperties = mediaPropertiesFunc(parameters)
  mediaProperties = MCmatlab.mediumProperties;

  j=1;
  mediaProperties(j).name  = 'Air';
  mediaProperties(j).mua   = 1e-8; % Absorption coefficient in cm^-1
  mediaProperties(j).mus   = 1e-8; % Scattering coefficient in cm^-1
  mediaProperties(j).g     = 1; % Henyey-Greenstein scattering anisotropy
  mediaProperties(j).n     = 1; % Refractive index

  j=2;
  lambda = 660; % wavelength in nm
  SC_Vw = 0.05; % fraction of water
  mua_w_660 = 0.0032; % absorption coefficient for water at 660 nm in cm^-1
  mua_baseline = 7.84*10^7*lambda^-3.255; 
  mediaProperties(j).name  = 'Stratum Corneum';
% mediaProperties(j).mua   = ((0.1-0.3*10^-4*lambda)+0.125*mua_baseline)*(1-SC_Vw)*mua_w_660; % Absorption coefficient in cm^-1
  mediaProperties(j).mua   = 2.5; % Absorption coefficient in cm^-1
  mediaProperties(j).mus   = 54.66; % Scattering coefficient in cm^-1
  mediaProperties(j).g     = 0.942; % Henyey-Greenstein scattering anisotropy
  mediaProperties(j).n     = 1.5; % Refractive index

  j=3;
  LE_Vmel = 0.11; % fraction of melanin
  LE_Vw = 0.2; % fraction of water
  mua_mel_660 = 269.4; % absorption coefficient for melanin at 660 nm in cm^-1
  mediaProperties(j).name  = 'Living Epidermis';
% mediaProperties(j).mua   = (LE_Vmel*mua_mel_660+(1-LE_Vmel)*mua_baseline)*(1-LE_Vw)+LE_Vw*mua_w_660; % Absorption coefficient in cm^-1
  mediaProperties(j).mua   = 2.5;
  mediaProperties(j).mus   = 54.66; % Scattering coefficient in cm^-1
  mediaProperties(j).g     = 0.8; % Henyey-Greenstein scattering anisotropy
  mediaProperties(j).n     = 1.34; % Refractive index

  j=4;
  V_A_SP = 0.17; % fraction of arterial blood
  V_V_SP = 0.17; % fraction of venous blood
  V_w_SP = 0.5; % fraction of water
  SaO2 = 0.97; % arterial oxygen saturation
  SvO2 = 0.67; % venous oxygen saturation
  mu_a_HbO2 = 1.4; % absorption coefficient oxygenated blood in cm^-1
  mu_a_HHb = 10; % absorption coefficient deoxygenated blood in cm^-1
  mu_a_A = SaO2*mu_a_HbO2+(1-SaO2)*mu_a_HHb;
  mu_a_V = SvO2*mu_a_HbO2+(1-SvO2)*mu_a_HHb;
  mediaProperties(j).name  = 'Stratum Papillare';
  mediaProperties(j).mua   = V_A_SP*mu_a_A+V_V_SP*mu_a_V+V_w_SP*mua_w_660+(1-V_A_SP+V_V_SP+V_w_SP)*mua_baseline; % Absorption coefficient in cm^-1
  mediaProperties(j).mus   = 53.62; % Scattering coefficient in cm^-1
  mediaProperties(j).g     = 0.9; % Henyey-Greenstein scattering anisotropy
  mediaProperties(j).n     = 1.4; % Refractive index

  j=5;
  V_A_SR = 0.07; % fraction of arterial blood
  V_V_SR = 0.07; % fraction of venous blood
  V_w_SR = 0.7; % fraction of water
  mediaProperties(j).name  = 'Stratum Reticulare';
  mediaProperties(j).mua   = V_A_SR*mu_a_A+V_V_SR*mu_a_V+V_w_SR*mua_w_660+(1-V_A_SR+V_V_SR+V_w_SR)*mua_baseline; % Absorption coefficient in cm^-1
  mediaProperties(j).mus   = 53.62; % Scattering coefficient in cm^-1
  mediaProperties(j).g     = 0.8; % Henyey-Greenstein scattering anisotropy
  mediaProperties(j).n     = 1.4; % Refractive index

  j=6;
  g = 0.75;
  mediaProperties(j).name  = 'Hypodermis';
  mediaProperties(j).mua   = 0.4; % Absorption coefficient in cm^-1
  mediaProperties(j).mus   = 12.1/(1-g); % Scattering coefficient in cm^-1
  mediaProperties(j).g     = 0.75; % Henyey-Greenstein scattering anisotropy
  mediaProperties(j).n     = 1.44; % Refractive index
 
  j=7;
  mediaProperties(j).name  = 'Muskel';
  mediaProperties(j).mua   = 1.28; % Absorption coefficient in cm^-1
  mediaProperties(j).mus   = 73.56; % Scattering coefficient in cm^-1
  mediaProperties(j).g     = 0.93; % Henyey-Greenstein scattering anisotropy
  mediaProperties(j).n     = 1.41; % Refractive index

  j=8;
  mediaProperties(j).name  = 'Knochen';
  mediaProperties(j).mua   = 0.667; % Absorption coefficient in cm^-1
  mediaProperties(j).mus   = 280; % Scattering coefficient in cm^-1
  mediaProperties(j).g     = 0.93; % Henyey-Greenstein scattering anisotropy
  mediaProperties(j).n     = 1.41; % Refractive index
end