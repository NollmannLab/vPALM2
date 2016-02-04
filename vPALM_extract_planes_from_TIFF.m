%% initialises
clear

window=[10, 10]; %window for extraction of planes.
piezo_step=50; % in nm

%% loads data
disp('Give the filename of the image file');
[fullFileName,path1] = uigetfile('*.tif','Open Tif File');
cd(path1);
[file.path,file.name, file.ext]=fileparts(fullFileName);
filename_output=strcat(file.name,'_roi.tif');
first_frame=1;
finfo = imfinfo(fullFileName);
Nframes=numel(finfo);
I = loadtiff(fullFileName);

%% defines ROI
info{1,1}.roi=[1 1 170 170 ]; 
info{1,2}.roi=[1 170 170 170 ]; 
info{1,3}.roi=[1 340 170 170 ]; 
info{2,1}.roi=[170 1 170 170 ]; 
info{2,2}.roi=[170 170 170 170 ]; 
info{2,3}.roi=[170 340 170 170 ]; 
info{3,1}.roi=[340 1 170 170 ]; 
info{3,2}.roi=[340 170 170 170 ]; 
info{3,3}.roi=[340 340 170 170 ]; 

for i=1:3
    for j=1:3
        info{i,j}.im_crop=zeros(info{i,j}.roi(3)+1,info{i,j}.roi(4)+1,Nframes);
        info{i,j}.intensity=[];
    end
end


%% makes ROI and calculates entropy
figure('Position',[800 200 1200 1000],'Color',[1 1 1],'name','All planes');
for iframe=first_frame:1:Nframes

    for i=1:3
        for j=1:3

            info{i,j}.im_crop(:,:,iframe)=imcrop(I(:,:,iframe),info{i,j}.roi); 
            info{i,j}.intensity=[info{i,j}.intensity, entropy(uint16(info{i,j}.im_crop(:,:,iframe)))];
            
        end
    end
    subplot(10,10,min([iframe 100]))
    imagesc(info{2,2}.im_crop(:,:,iframe));
    title(strcat('Frame #',num2str(iframe)))

iframe
end

%% analyses entropy profiles and saves
figure, hold on
colors=jet(9);

iplane=1;plane_position=[];

for i=1:3
    for j=1:3
        x=[first_frame:Nframes];
        plot(x*piezo_step,info{i,j}.intensity,'color',colors(iplane,:))
        
        [~,info{i,j}.maxfr]=min(info{i,j}.intensity);
        frange=find(x>max([info{i,j}.maxfr-2*window(1),1]) & x<min([Nframes,info{i,j}.maxfr+2*window(2)]));
        % refines center location by fitting a second order polynomial
        p=polyfit(x(frange),info{i,j}.intensity(frange),2);
        p2=[2*p(1) p(2)];
        if sum(p2<Inf)==2
            x0x2=roots(p2);
        else
            x0x2=0.0;
        end
        info{i,j}.focal_plane=(x0x2);
        line([ piezo_step*info{i,j}.focal_plane  piezo_step*info{i,j}.focal_plane],[min(info{i,j}.intensity)-1 max(info{i,j}.intensity)+1],'color',colors(iplane,:));
        disp(strcat('Focal plane [',num2str(i),num2str(j),'] #',num2str(iplane),': ',num2str(info{i,j}.focal_plane)))
        plane_position=[plane_position,info{i,j}.focal_plane];
        
        filename_new=strcat(file.name,'_roi',num2str(i),num2str(j),'.tif');
        saveastiff(uint16(info{i,j}.im_crop(:,:,max([1 info{i,j}.maxfr-window(1)]):min([Nframes info{i,j}.maxfr+window(2)]))),filename_new );
        iplane=iplane+1;
    end
end
plane_positions=(sort(plane_position))*piezo_step

plane_distances=diff(sort(plane_position))*piezo_step
disp(strcat('Mean plane distance=',num2str(mean(plane_distances)),'+-',num2str(std(plane_distances))));

title('Entropy versus pzt position for 9 planes')
xlabel('pzt position, nm');
ylabel('Entropy, a.u.');

