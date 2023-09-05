function plot_pd_wrinkels()
    
    angle = [0,5,10,15,20,25,30,35,40,45];
    
    figure
    for depth = 0:0.02:0.1
        pd = zeros(length(angle),1);
        i = 1;
        for a = angle
            pd(i) = get_fluence_wrinkels(a,depth);
            pd(i) = pd(i)*(0.53/1000);
            i = i+1;
        end
        hold on
        plot(angle,pd,'DisplayName',num2str(depth));
        xlabel('polar angle of the source (deg)')
        ylabel('penetration depth');
        set(gca,'YDir','reverse');
        legend;
        title('Penetration depth of decected ligh for ippg with changing angle');
    end
    hold off

end



function PD = get_fluence_wrinkels(angle,depth)
    PLOTON = true;
    % Laden der Daten
    filename = sprintf('outputs/ippg/%ddeg470um2_%gwrinkels.mat',angle,depth);
    fluence = load(filename,'fluence');
    fluence = fluence.fluence;
    geometry = load(filename,'geometry');
    geometry = geometry.geometry;
    fluence(geometry == 1) = 0;
    
    % Berechnung Flux(2D)
    flux = zeros([size(fluence,1),size(fluence,3)]);
    for i = 1:size(fluence,3)
        for j = 1:size(fluence,1)
            flux(j,i) = sum(fluence(j,:,i));
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
