function im = roughness(type,intersection,dimxy,mean)
   
    if (type == 1)%random
        
        filtersize = dimxy/10;
        X = ones(dimxy+filtersize-1,dimxy+filtersize-1);
        H = ones(filtersize,filtersize)./filtersize^2;
        if(intersection == 1)%air-sc
        %mean = 20;
        std = 3;
        elseif(intersection == 2)%sc-le
        %mean = 40;
        std = 4;
        elseif(intersection == 3)%le-pd
        %mean = 140;
        std = 3;
        end
        r = normrnd(0,std,size(X));
        im = filter2(H,r,'valid');
        s = std2(im);
        im = im.*(std/s);
        im = im+mean;
        im = im.*10e-5;
        figure 
        mesh(im);

    elseif (type == 2)%sulci cutis
        
        if(intersection == 1)%air-sc
        %mean = 20;
        ra = 14.2;
        rz = 75.7;

        std = 3;
        im = zeros(dimxy,dimxy);
        figure
        while(1)
        winkel = randsample(180,1);
        transly = randsample(400,1);
        x = linspace(1,200,200);
        y = tan((winkel*2*pi/360))*x+transly;
        y = cast(y,'int16');
        y(y<1) = 1;
        y(y>200) = 200;
        line = abs(normrnd(ra,rz,[1,dimxy]));
        line(y==1) = 0;
        line(y==200) = 0;
        line(line >= rz) = rz;
        for i = 1:200
            if(im(x(i),y(i)) == 0)
                im(x(i),y(i)) = im(x(i),y(i)) + line(i);
            else
                im(x(i),y(i)) = 0;
            end
        end
        im2 = imgaussfilt(im,std,"FilterSize",[7,7]);
        mesh(im2);
        ra_i = 1/(dimxy*dimxy) * sum(sum(im2));
        if (ra_i > ra)
            break;
        end
        end
        RoughnessStore.setgetIm(im);
        im = imgaussfilt(im,std,"FilterSize",[7,7]);

        elseif(intersection == 2)%sc-le
        %mean = 40;
        ra = 14.2;
        rz = 75.7;
        std = 3;
        im = RoughnessStore.setgetIm;
        im = imgaussfilt(im,std,"FilterSize",[7,7]);

        elseif(intersection == 3)%le-pd
        %mean = 140;
        ra = 15.3;
        rz = 81.4;
        std = 3;
        im = RoughnessStore.setgetIm;
        im = im.*(15.3/14.2); 
        im = imgaussfilt(im,std,"FilterSize",[9,9]);
        end
        
        %s = std2(im);
        im = im+mean;
        %im([1, 200],:) = 0;
        %im(:,[1, 200]) = 0;
        im = im.*10e-5;
        figure
        mesh(im);
    
    
    elseif (type == 3)%sinus
        if(intersection == 1)%air-sc
            A = 2;
            omegax = 2*pi*1/(100/2);
            omegay = 2*pi*1/(150/2);
            %mean = 20;
        elseif(intersection == 2)%sc-le
            A = 2.5;
            omegax = 2*pi*1/(80/2);
            omegay = 2*pi*1/(80/2);
            %mean = 40;
        elseif(intersection == 3)%le-pd
            A = 20;
            omegax = 2*pi*1/(50/2);
            omegay = 2*pi*1/(45/2);
            %mean = 140;
        end
        
        offsetx = randsample(20,1);
        offsety = randsample(20,1);
        figure
        im1 = zeros(dimxy,dimxy);
        im2 = zeros(dimxy,dimxy);
        for x = 1:dimxy
            x_i = x*0.1*10e3/dimxy;
            im1(x,:) = A*sin(x_i*omegax+offsetx); 
        end
        for y = 1:dimxy
            y_i = y*0.1*10e3/dimxy;
            im2(:,y) = (A*sin(y_i*omegay+offsety));
        end
        im = im1.*im2;
        im = im+mean;
        im = im.*10e-5;
        mesh(im);
    end
end