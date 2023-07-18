function PD = plot_fluence()
    fluence = load('outputs\ppg\15deg470um4.mat','fluence');
    fluence = fluence.fluence;
    flux = zeros([size(fluence,2),size(fluence,3)]);
    for i = 1:size(fluence,3)
        for j = 1:size(fluence,2)
            flux(j,i) = sum(fluence(:,j,i));
        end
    end
    flux = flux./max(max(flux));
    flux_plt = rot90(rot90(rot90((flux))));
    figure
    colormap("hot")
    imshow(flux_plt);
    
    
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
