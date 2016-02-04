function [PALM_image] = vPALM2_analysis_image_reconstruction_v4b(X, Y, s, A, ps, border, Pixel_to_photon,uniform_peaks, max_intensity, hot, Contrast_tool)

% Keep in mind that X and Y are definied as the two axis of a cartesian
% referencial. For the image, you need to switch from this cartesian
% notation to a matricial one using rows and columns. In that case, remember that X is related to 
% columns and Height-Y to the rows.


N = size(X,1);

overflowed = find(A>max_intensity); % Look for all the events that display an intensity greater than max_intensity
if uniform_peaks
    A(overflowed) = max_intensity; % Allow us to control the contrast of the image
end

Xmin = min(X);
Ymin = min(Y);
Xmax = max(X);
Ymax = max(Y);

% The size of the pixel is defined by 'ps' and its default value is set at 5nm=5e-3?m. 
 
xx1 = floor(Xmin/ps)*ps:ps:ceil(Xmax/ps)*ps;
xxcenters = xx1(1)+ps/2:ps:xx1(end)-ps/2;
yy1 = floor(Ymin/ps)*ps:ps:ceil(Ymax/ps)*ps;
yycenters = yy1(1)+ps/2:ps:yy1(end)-ps/2;

imsize_half=5;
imsize=imsize_half*2+1;
        
Nx = round (((Xmax+border) - (Xmin-border))/ps); %Nx and Ny must be integers 
Ny = round(((Ymax+border) - (Ymin-border))/ps);

hwb = waitbar(0,['Calculating image with ',num2str(size(xx1,2)*size(yy1,2)),' points...']);

I = zeros (Ny+10,Nx+10); %creates a matrix Ny x Nx filled with zeros - Note that x is associated to the columns and y to the rows
        Imax=0;     
        for i=1 : N
            waitbar(i/N);
            
            [Ispot, Nxspot] = gaussian_spot_vPALM(s(i),ps);
            
%             [Ispot, Nxspot] = gaussian_spot_vPALM(s(i),imsize)
%             size_spot=size(Ispot);
        
            [aux,i1] = min(abs(yycenters-Y(i)));
            [aux,j1] = min(abs((xxcenters-X(i))));
            iL = i1-(Nxspot-1)/2;
            jL = j1-(Nxspot-1)/2;
            iR = i1 + (Nxspot-1)/2;
            jR = j1 + (Nxspot-1)/2;
            iL2 = max([iL 1]) ;
            jL2 = max([jL 1]) ;
            iR2 = min([iR Ny]) ;
            jR2 = min([jR Nx]) ;
            % set the intensity of the spot
            if uniform_peaks
                spot_intensity = max_intensity;
            else
                spot_intensity = A(i)*Pixel_to_photon;
            end
            % add the gaussian spot to the image
            I(iL2:iR2, jL2:jR2) = I(iL2:iR2, jL2:jR2)  + Ispot(1+iL2-iL:end-(iR-iR2), 1+jL2-jL:end-(jR-jR2)) * spot_intensity  ;
            tmpImax=max(max(I(iL2:iR2, jL2:jR2)));
            if Imax<tmpImax
                Imax=tmpImax;
            end
         end


Imin = min(I(:));

h1 = figure;
hold on

N = 1000/ps;
for k = 1 : N
%     I(10+k,10:10+round(N/10)) = Imax; %Scale bar added here (1?m)
    I(10:10+round(N/10),10+k) = Imax; %Scale bar added here (1?m)

end

yy1 = yy1 - min(yy1);
xx1 = xx1 - min(xx1);

if hot
    I_normalized0 = round(16384/(Imax-Imin)*(I - Imax) + 16384);
%     I_normalized = flipud(I_normalized);
% gaussian filter
    h = fspecial('gaussian', [3 3], 2);
    I_normalized=imfilter(I_normalized0,h);
%     imagesc(yy1,xx1,I_normalized);
    imagesc(xx1,yy1,(I_normalized));

    xlabel('X (nm)');
    ylabel('Y (nm)');  
    axis equal
    colormap('hot');
    colorbar();
    if Contrast_tool
        imcontrast(h1);
    end
else
%     imagesc(yy1,xx1,I_normalized); 
    imagesc(xx1,yy1,fliplr((I_normalized')));

    axis equal
    xlabel('X (nm)');
    ylabel('Y (nm)');  
    if Contrast_tool
        imcontrast(h1);
    end
end

% PALM_image = I_normalized';
PALM_image = I_normalized;

waitbar();

Profile = questdlg('Do you want to calculate a line profile?','Line Profile','Yes','No','Yes');
switch Profile
    case 'Yes'
        
        Proceed = 1;
        
        while Proceed
            
            h = imline(gca(h1));
            position = wait(h); % wait until you double-click on the image and return the position of the line
            
            xline = position(:,1);
            x2 = xline(2);
            x1 = xline(1);
            yline = position(:,2);
            y2 = yline(2);
            y1 = yline(1);
            
            npoints = round(sqrt((xline(2)-xline(1))^2+(yline(2)-yline(1))^2));
            
            x = [xline(1):(xline(2)-xline(1))/npoints : xline(2)];
            y = [yline(1):(yline(2)-yline(1))/npoints : yline(2)];
            
            hold on,
            
            Iline = zeros(npoints,1);
            Iline2 = zeros(npoints,1);
            np = 2;
            
            for i = 1 : npoints+1
%                 i
%                 y(i)/ps
%                 x(i)/ps
                Iline(i) = mean(mean(I_normalized(round(y(i)/ps) - np : round(y(i)/ps) + np, round(x(i)/ps) - np : round(x(i)/ps) + np)));
                Iline2(i) = I_normalized(round(y(i)/ps),round(x(i)/ps));
            end
            distances = sqrt((x-x1).^2+(y-y1).^2);
            
            figure
            plot(distances,Iline2,'-g')
            hold on
            plot(distances,Iline,'-r')
            
            Proceed = questdlg('Do you want to calculate another profile?','Another','Yes','No','Yes');
            switch Proceed
                case 'Yes'
                    Proceed = 1;
                case 'No'
                    Proceed = 0;
            end
        end
end
;