%% sliders setup
function slider_setup(h,m,parameters)

% if h.parameters.fit_3dxcorr==1
%     set(h.Imin_slider,'Max',max(m(:,12)),'Min', min(m(:,12)),'Value',min(m(:,12)));
%     set(h.Imax_slider,'Max',max(m(:,12)),'Min', min(m(:,12)),'Value',max(m(:,12)));
%     set(h.Imax,'String',num2str(max(m(:,12))));
%     set(h.Imin,'String',num2str(min(m(:,12))));    
% else
set(h.h_sliders.Imin_slider,'Max',max(m(:,4)),'Min', min(m(:,4)),'Value',min(m(:,4)));
set(h.h_sliders.Imax_slider,'Max',max(m(:,4)),'Min', min(m(:,4)),'Value',max(m(:,4)));
set(h.h_sliders.Imax,'String',num2str(max(m(:,4))));
set(h.h_sliders.Imin,'String',num2str(min(m(:,4))));
% end

set(h.h_sliders.Frmin_slider,'Max',max(m(:,1)),'Min', min(m(:,1)),'Value',min(m(:,1)));
set(h.h_sliders.Frmax_slider,'Max',max(m(:,1)),'Min', min(m(:,1)),'Value',max(m(:,1)));
set(h.h_sliders.Frmax,'String',num2str(max(m(:,1))));
set(h.h_sliders.Frmin,'String',num2str(min(m(:,1))));

set(h.h_sliders.Xmin_slider,'Max',max(m(:,2)),'Min', min(m(:,2)),'Value',min(m(:,2)));
set(h.h_sliders.Xmax_slider,'Max',max(m(:,2)),'Min', min(m(:,2)),'Value',max(m(:,2)));
set(h.h_sliders.Xmax,'String',num2str(max(m(:,2))));
set(h.h_sliders.Xmin,'String',num2str(min(m(:,2))));

set(h.h_sliders.Ymin_slider,'Max',max(m(:,3)),'Min', min(m(:,3)),'Value',min(m(:,3)));
set(h.h_sliders.Ymax_slider,'Max',max(m(:,3)),'Min', min(m(:,3)),'Value',max(m(:,3)));
set(h.h_sliders.Ymax,'String',num2str(max(m(:,3))));
set(h.h_sliders.Ymin,'String',num2str(min(m(:,3))));


if parameters.fit_3dxcorr==1 %      plots 3d xcorr localizations

    set(h.h_sliders.Zmin_slider,'Max',max(m(:,13)),'Min', min(m(:,13)),'Value',min(m(:,13)));
    set(h.h_sliders.Zmax_slider,'Max',max(m(:,13)),'Min', min(m(:,13)),'Value',max(m(:,13)));
    set(h.h_sliders.Zmax,'String',num2str(max(m(:,13))));
    set(h.h_sliders.Zmin,'String',num2str(min(m(:,13))));

elseif parameters.fit_3dgaussian==1 %      plots 3d gaussian localizations
%     set(h.Zmin_slider,'Max',max(m(:,14)./m(:,13)),'Min', min(m(:,14)./m(:,13)),'Value',min(m(:,14)./m(:,13)));
% %     set(h.Zmax_slider,'Max',max(m(:,14)./m(:,13)),'Min', min(m(:,14)./m(:,13)),'Value',max(m(:,14)./m(:,13)));
    set(h.h_sliders.Zmax_slider,'Max',10000,'Min', -10000,'Value',0);
    set(h.h_sliders.Zmin_slider,'Max',10000,'Min', -10000,'Value',0);

%     set(h.h_sliders.Zmax,'String',num2str(max(m(:,14)./m(:,13))));
    set(h.h_sliders.Zmax,'String',num2str(10000));

%     set(h.Zmin,'String',num2str(min(m(:,14)./m(:,13))));
    set(h.h_sliders.Zmin,'String',num2str(-10000));

    
else
    ;
    
end