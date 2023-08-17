function plot_pd(roughness_type,sensor_type)
    switch roughness_type
        case 'random'
            r = 1;
        case 'sulci cutis'
            r = 2;
        case 'sinus'
            r = 3;
        case 'normal'
            r = 4;
        otherwise
            r = 1;
    end
    nmLib = load('spectralLIBskin.mat','nmLIB');
    nmLib = nmLib.nmLIB;
    angle = [0,5,10,15,20,25,30,35,40,45];
    
    figure
    for nm = nmLib
        pd = zeros(length(angle),1);
        i = 1;
        for a = angle
            pd(i) = plot_fluence(sensor_type,a,nm,r);
            pd(i) = pd(i)*(0.53/1000);
            i = i+1;
        end
        hold on
        plot(angle,pd,'DisplayName',num2str(nm));
        xlabel('polar angle of the source (deg)')
        ylabel('penetration depth');
        set(gca,'YDir','reverse');
        legend;
        title(sprintf('Penetration depth of decected ligh for %s with changing angle',sensor_type));
    end
    hold off

end
