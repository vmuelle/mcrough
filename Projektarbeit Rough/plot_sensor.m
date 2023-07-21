function plot_sensor(roughness_type,sensor_type)
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
    
    for nm = nmLib
        intensity = zeros(length(angle),1);
        i = 1;
        for a = angle
            filename = sprintf('outputs/%s/%ddeg%dum%d.mat',sensor_type,a,nm,r);
            sensor_image = load(filename,'image');
            sensor_image = sensor_image.image;
            %figure
            %imshow(sensor_image);
            intensity(i) = sum(sum(sensor_image));
            i = i+1;
        end
        figure
        plot(angle,intensity);
        xlabel('polar angle of the source (deg)')
        ylabel('detected energy');
        title(sprintf('Detected energy for %s with changing angle with light at %d nm',sensor_type,nm));
    end
end