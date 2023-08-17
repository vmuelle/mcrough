function PD = plot_fluence(sensor_type,angle,nm,roughness_type)

    PLOTON = false;
    
    % Laden der Daten
    fluence = load(sprintf('outputs/%s/%ddeg%dum%d.mat',sensor_type,angle,nm,roughness_type),'fluence');
    fluence = fluence.fluence;
    geometry = load(sprintf('outputs/%s/%ddeg%dum%d.mat',sensor_type,angle,nm,roughness_type),'geometry');
    geometry = geometry.geometry;
    fluence(geometry == 1) = 0;
    
    % Berechnung Flux(2D)
    flux = zeros([size(fluence,2),size(fluence,3)]);
    for i = 1:size(fluence,3)
        for j = 1:size(fluence,2)
            flux(j,i) = sum(fluence(:,j,i));
        end
    end
    flux = flux./max(max(flux));
    if(PLOTON)
        flux_plt = rot90(rot90(rot90((flux))));
        figure
        colormap("hot")
        imshow(flux_plt);
    end

    % Brechnung Fluence(1D)
    z = linspace(0,0.53,1000)';
    fl = zeros(length(z),1);
    for i=1:length(z)
        fl(i,1) = sum(flux(:,i),'all');
    end
    fluence_norm = fl./max(fl);
    if(PLOTON)
        % Plot
        figure
        plot(z*10,fluence_norm);
        xlabel('Hauttiefe in mm');
        ylabel('a.u.');
        grid minor;
    end
    
    % Berechnung Penetration Depth
    sum_flux = sum(sum(flux));
    sum_flux_i = 0;
    for i = 1:size(flux,2)
        sum_flux_i = sum_flux_i+sum(flux(:,i)); 
        if(sum_flux_i >= sum_flux * 0.632)
            PD = i;
            break;
        end
    end
end
