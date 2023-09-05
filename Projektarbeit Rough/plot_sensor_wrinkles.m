function plot_sensor_wrinkels()
   
  
    angle = [0,5,10,15,20,25,30,35,40,45];
    
    figure
    for depth = 0:0.02:0.1
        intensity = zeros(length(angle),1);
        i = 1;
        for a = angle
            filename = sprintf('outputs/ippg/%ddeg470um2_%gwrinkels.mat',a,depth);
            sensor_image = load(filename,'image');
            sensor_image = sensor_image.image;
            %figure
            %imshow(sensor_image);
            intensity(i) = sum(sum(sensor_image));
            i = i+1;
        end
        hold on
        plot(angle,intensity,'DisplayName',num2str(depth));
        %set(gca,'YDir','reverse');
        legend;
        xlabel('polar angle of the source (deg)')
        ylabel('detected energy');
        title('Detected energy for ippg with changing angle');
    end
    hold off
end