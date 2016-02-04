function vPALM_plot_reconstruction(hObject, eventdata, h)
status(h,'g','Calculating reconstruction...');
m=h.m;
% pointing_precision_px=str2num(get(h.h_sliders.detect_thresh,'String'))./str2num(get(h.pixelsize,'String'));

% pointing_precision_px=str2num(get(h.detect_thresh,'String'))./str2num(get(h.pixelsize,'String'));
pointing_precision_px=str2num(get(h.resolution_reconstruction,'String'))./str2num(get(h.pixelsize,'String'));
ps=.05; %pixel size in px units
logyes=0 ;% 1 for log

m=filter_sliders(h,h.m);
x=m(:,2);
y=m(:,3);
A=m(:,4);
t=m(:,1);
if get(h.chb_fix_resolution,'Value')==1
    if h.parameters.fit_2dgaussian==1 | (h.parameters.fit_3dgaussian==0 & h.parameters.fit_3dxcorr==0)
    %     s=m(:,5);
        s=ones(size(m(:,2),1),1)*pointing_precision_px*str2num(get(h.pixelsize,'String'));
    elseif h.parameters.fit_3dgaussian==1 
    %     s=sqrt(m(:,13)./m(:,14));
        s=ones(size(m(:,2),1),1)*pointing_precision_px*str2num(get(h.pixelsize,'String'));
    elseif h.parameters.fit_3dxcorr==1
        s=ones(size(m(:,2),1),1)*pointing_precision_px*str2num(get(h.pixelsize,'String'));
    else
    %     s=sqrt(m(:,5));
        s=ones(size(m(:,2),1),1)*pointing_precision_px*str2num(get(h.pixelsize,'String'));

    end
else
    switch get(h.load_options,'Value')
    
    case 1 %qPALM

    case 2 % PALMcbs
        
    
    case 4 %MTT
      
        s=ones(size(m(:,2),1),1)*pointing_precision_px*str2num(get(h.pixelsize,'String'));

    case 5 %ENS
       
    case 6 %micromanager
        s=m(:,5);%./str2num(get(h.pixelsize,'String'));

    case 7 %micromanager
        
        s=m(:,5);%./str2num(get(h.pixelsize,'String'));
   case 8 %micromanager
        
        s=m(:,5);%./str2num(get(h.pixelsize,'String'));

        
    end
end

if get(h.chb_DIPimage,'Value')==0
%     px = str2num(get(h.pixelsize,'String'));
%     x = x*px;
%     y = y*px;
%     uniform_peaks=get(h.chb_uniform_peaks,'Value');
%     border=10;
%     Pixel_to_photon=1;
%     max_intensity=1;
%     hot=1;
%     Contrast_tool=1;
%     ps=str2num(get(h.stepsize_reconstruction,'String'));
%     % [I,xcoor,ycoor,Imax]=vPALM_reconstruction_v5(x,y,s,A, ps,10,1,1,1,1,logyes);
%     [PALM_image_no_gaussian]=vPALM2_analysis_image_reconstruction_v4b(x, y, s, A, ps, border, Pixel_to_photon,uniform_peaks, max_intensity, hot, Contrast_tool);
disp('This functionality is no longer supported')
else
    % DIPimage
%     ps=str2num(get(h.stepsize_reconstruction,'String'));
%     px = str2num(get(h.pixelsize,'String'));
%     
%     coords(:,1) = x-min(x);
%     coords(:,2) = y-min(y);
%     coords(:,3) = t;
%     loc(:,1) = x-min(x);
%     loc(:,2) = y-min(y);
% 
% %     uniform_peaks=get(h.chb_uniform_peaks,'Value');
%     s= pointing_precision_px*str2num(get(h.pixelsize,'String'));
% 
% %     
% %     
%     superzoom = ceil(px/ps);
%     szx = superzoom * ceil(max(coords(:,1)));
%     szy = superzoom * ceil(max(coords(:,2)));
% % 
%     im = binlocalizations(coords, szx, szy, superzoom);
%     h1=dipshow(im);
%     dipmapping(h1,[0 5]);%,'colormap',hot);
%% 2d histogram + gaussian bluring... easy and fast
ps=str2num(get(h.stepsize_reconstruction,'String'));
px = str2num(get(h.pixelsize,'String'));
s= pointing_precision_px*str2num(get(h.pixelsize,'String'));

vPALM2_reconstruction_time(x,y,t,s,ps,px);

end
% axes(h.axes2)
% Imin=0;
% % f1=figure('Color',[0 0 0 ]);
% % I = round(16384/(Imax-Imin)*(I - Imax) + 16384);
% imagesc(I)
% %     if logyes
% %         [C,h]=contour(flipud(rot90(log((I+100)))));
% %     else
% %         [C,h]=contour(flipud(rot90(((I+100)))));
% %     end
%     xlabel('X (\mum)','color',[.8 .8 .8]);
%     ylabel('Y (\mum)','color',[.8 .8 .8]);  
% %     axis equal
% %     colormap('hot');
%     set(gca,'Xcolor',[0.5 0.5 0.5]);
%     set(gca,'Ycolor',[0.5 0.5 0.5]);
%     set(gca,'Color',[0 0 0]);

status(h,'g','done');
