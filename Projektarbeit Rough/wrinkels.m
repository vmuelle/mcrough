close all
spLIB = load('spectralLIBskin.mat');
model1 = MCmatlab.model;
model1.G.nx                = 400; % Number of bins in the x direction
model1.G.ny                = 400; % Number of bins in the y direction
%model1.G.nz                = 400; % Number of bins in the z direction
model1.G.nz                = 1000;
model1.G.Lx                = 0.5; % x size of simulation cuboid in cm
model1.G.Ly                = 0.5; % y size of simulation cuboid in cm
%model1.G.Lz                = 0.037; % z size of simulation cuboid in cm
%model1.G.Lz                = 0.822;
model1.G.Lz                = 0.53;

%% Monte Carlo simulation
model1.MC.simulationTimeRequested  = 5; % Time duration of the simulation in min
model1.MC.nExamplePaths            = 100;
model1.MC.matchedInterfaces        = false;
model1.MC.smoothingLengthScale     = model1.G.Lx*10;
model1.MC.boundaryType             = 1; % 1: All cuboid boundaries are escaping
model1.MC.wavelength               = 470; % Excitation wavelength in nm, used for determination of optical properties for excitation light

roughness_type = 2;
nm = 1;

model1.MC.wavelength         = spLIB.nmLIB(nm); % Excitation wavelength in nm, used for determination of optical properties for excitation light

        
% 
% model1.G.geomFuncParams = {roughness_type,0.04};
% model1.G.geomFunc = @geometryDefinition;
% model1.G.mediaPropParams = {spLIB,nm};
% model1.G.mediaPropertiesFunc = @mediaPropertiesFunc; % Media properties defined as a function at the end of this file
% model1 = plot(model1,'G');
% ippg_wrinkels(45,model1,roughness_type,0.04);


for winkel = 0:5:45
    for depth = 0:0.02:0.1
        model1.G.geomFuncParams = {roughness_type,depth};
        model1.G.geomFunc = @geometryDefinition;
        model1.G.mediaPropParams = {spLIB,nm};
        model1.G.mediaPropertiesFunc = @mediaPropertiesFunc; % Media properties defined as a function at the end of this file
        model1 = plot(model1,'G');
        ippg_wrinkels(winkel,model1,roughness_type,depth);
    end
end

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
    %zsurf = 20;
    %SC_thick = 40;
    %LE_thick = 140;
    PD_thick = 0.015;
    RD_thick = 0.15;
    Hypo_thick = 0.1;
    Muskel_thick = 0.225;
    wrinkel_depth_val = parameters{2};

    dimxy = size(X,1); 
    
%     im1 = roughness(parameters{1},1,dimxy,zsurf*1e4);
%     im2 = roughness(parameters{1},2,dimxy,zsurf*1e4+SC_thick*1e4);
%     im3 = roughness(parameters{1},3,dimxy,zsurf*1e4+SC_thick*1e4+LE_thick*1e4);
%     [X,Y] = meshgrid(0:0.5/400:0.5-(0.5/400),0:0.5/400:0.5-(0.5/400));
%     X = X-0.25;
%     Y = Y-0.25;
% 
%     figure
%     mesh(X,Y,im1,'FaceColor','red','DisplayName','Oberfläche','EdgeColor',[0.7,0.7,0.7]);
%     hold on
%     mesh(X,Y,im2,'FaceColor','blue','DisplayName','SC-LE Übergang','EdgeColor',[0.7,0.7,0.7]);
%     mesh(X,Y,im3,'FaceColor','green','DisplayName','LE-PD Übergang','EdgeColor',[0.7,0.7,0.7]);
%     xlabel('X(cm)')
%     ylabel('Y(cm)')
%     zlabel('Z(cm)')
%     legend()
%     set(gca,'ZDir','reverse')
%     hold off

    M = ones(size(X)); % air
    im1 = wrinkel_depth(wrinkel_depth_val,dimxy);
    im2 = roughness(parameters{1},1,dimxy,zsurf*1e4);
    im2 = im2+im1;
    M(Z > im2) = 2; % stratum corneum
    im3 = roughness(parameters{1},2,dimxy,zsurf*1e4+SC_thick*1e4);
    im3 = im3+im1;
    M(Z > im3) = 3;
    im4 = roughness(parameters{1},3,dimxy,zsurf*1e4+SC_thick*1e4+LE_thick*1e4);
    im4 = im4+im1;
    M(Z > im4) = 4;
    im5 = ones(dimxy,dimxy).*(zsurf+SC_thick+LE_thick+PD_thick);
    im5 = im5+im1;
    M(Z > im5) = 5;
    im6 = ones(dimxy,dimxy).*(zsurf+SC_thick+LE_thick+PD_thick+RD_thick);
    im6 = im6+im1;
    M(Z > im6) = 6;
    im7 = ones(dimxy,dimxy).*(zsurf+SC_thick+LE_thick+PD_thick+RD_thick+Hypo_thick);
    im7 = im7+im1;
    M(Z > im7) = 7;
    M(Z > (zsurf+SC_thick+LE_thick+PD_thick+RD_thick+Hypo_thick+Muskel_thick)) = 8;
    %M(Z < wrinkel_depth(wrinkel_depth_val,dimxy,zsurf)) = 1;
    figure
    mesh(im6);
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

  lib = parameters{1};
  i = parameters{2};

  j=1;
  mediaProperties(j).name  = 'Air';
  mediaProperties(j).mua   = 1e-8; % Absorption coefficient in cm^-1
  mediaProperties(j).mus   = 1e-8; % Scattering coefficient in cm^-1
  mediaProperties(j).g     = 1; % Henyey-Greenstein scattering anisotropy
  mediaProperties(j).n     = 1; % Refractive index

  j=2;
  mediaProperties(j).name  = 'Stratum Corneum';
  mediaProperties(j).mua   = lib.mua_sc(i); % Absorption coefficient in cm^-1
  mediaProperties(j).mus   = lib.mus_sc(i);  % Scattering coefficient in cm^-1
  mediaProperties(j).g     = lib.g_sc(i);  % Henyey-Greenstein scattering anisotropy
  mediaProperties(j).n     = lib.n_sc(i);  % Refractive index

  j=3;
  mediaProperties(j).name  = 'Living Epidermis';
  mediaProperties(j).mua   = lib.mua_le(i);  % Absorption coefficient in cm^-1
  mediaProperties(j).mus   = lib.mus_le(i);  % Scattering coefficient in cm^-1
  mediaProperties(j).g     = lib.g_le(i);  % Henyey-Greenstein scattering anisotropy
  mediaProperties(j).n     = lib.n_le(i);  % Refractive index

  j=4;
  mediaProperties(j).name  = 'Stratum Papillare';
  mediaProperties(j).mua   = lib.mua_sp(i); % Absorption coefficient in cm^-1
  mediaProperties(j).mus   = lib.mus_sp(i);  % Scattering coefficient in cm^-1
  mediaProperties(j).g     = lib.g_sp(i);  % Henyey-Greenstein scattering anisotropy
  mediaProperties(j).n     = lib.n_sp(i);  % Refractive index

  j=5;
  mediaProperties(j).name  = 'Stratum Reticulare';
  mediaProperties(j).mua   = lib.mua_sr(i); % Absorption coefficient in cm^-1
  mediaProperties(j).mus   = lib.mus_sr(i); % Scattering coefficient in cm^-1
  mediaProperties(j).g     = lib.g_sr(i); % Henyey-Greenstein scattering anisotropy
  mediaProperties(j).n     = lib.n_sr(i); % Refractive index

  j=6;
  mediaProperties(j).name  = 'Hypodermis';
  mediaProperties(j).mua   = lib.mua_hd(i); % Absorption coefficient in cm^-1
  mediaProperties(j).mus   = lib.mus_hd(i); % Scattering coefficient in cm^-1
  mediaProperties(j).g     = lib.g_hd(i); % Henyey-Greenstein scattering anisotropy
  mediaProperties(j).n     = lib.n_hd(i); % Refractive index

  j=7;
  mediaProperties(j).name  = 'Muskel';
  mediaProperties(j).mua   = lib.mua_muskel(i); % Absorption coefficient in cm^-1
  mediaProperties(j).mus   = lib.mus_muskel(i); % Scattering coefficient in cm^-1
  mediaProperties(j).g     = lib.g_muskel(i); % Henyey-Greenstein scattering anisotropy
  mediaProperties(j).n     = lib.n_muskel(i); % Refractive index

  j=8;
  mediaProperties(j).name  = 'Knochen';
  mediaProperties(j).mua   = lib.mua_knochen(i); % Absorption coefficient in cm^-1
  mediaProperties(j).mus   = lib.mus_knochen(i);  % Scattering coefficient in cm^-1
  mediaProperties(j).g     = lib.g_knochen(i);  % Henyey-Greenstein scattering anisotropy
  mediaProperties(j).n     = lib.n_knochen(i);  % Refractive index

end

function im = wrinkel_depth(depth,dimxy)

    im = zeros(dimxy,dimxy);
    transly = dimxy/2;
    y = linspace(1,dimxy,dimxy);
    x = ones(dimxy).*transly;
    
    im(x,y) = 50;
    std = 20;
    im = imgaussfilt(im,std,"FilterSize",[101,101]);
    im = im.*depth;

    figure
    mesh(im);
  
end

function ippg_wrinkels(winkel, model,roughness_type,depth)

model.MC.lightSource.sourceType   = 2; % 5: X/Y factorizable beam (e.g., a rectangular LED emitter)

model.MC.lightSource.theta     = (winkel*2*pi)/360; % [rad] Polar angle of direction the light collector is facing
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

model.MC.depositionCriteria.onlyCollected = true;

model = runMonteCarlo(model);
model = plot(model,'MC');
fluence = model.MC.normalizedFluenceRate(:,:,:);
transmittance = model.MC.normalizedIrradiance_zpos;
reflectance = model.MC.normalizedIrradiance_zneg;
image = model.MC.lightCollector.image;
geometry = model.G.M_raw;
name = strcat('outputs/ippg/',string(winkel),'deg',string(model.MC.wavelength),'um',string(roughness_type),'_',string(depth),'wrinkels.mat');
save(name,'fluence','transmittance','reflectance','image','geometry');
%MCmatlab.plotAzFz(model1);

end