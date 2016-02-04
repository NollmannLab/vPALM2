%% plot_localizations_2d
function h=vPALM_plot_localizations_2d(hObject,eventdata,h)
tic
status(h,'r','working...');

if get(h.chb_loc_by_moment,'Value')==1
    xcol=2;
    ycol=3;
else
    xcol=10;
    ycol=11;
end


if get(h.chb_data_loaded,'Value')==1
    
    m0=h.m;
    parameters=h.parameters;
    v=get(h.Zcolormap,'Value');
       switch v
        case 1
            parameters.colormap='jet';
            parameters.colormap_label='resolution';
        case 2
            parameters.colormap='hot';
            parameters.colormap_label='z, nm';
        case 3
            parameters.colormap='hsv';
            parameters.colormap_label='density';
        case 4
            parameters.colormap='gray';
            parameters.colormap_label='ellipticity/t-density';
        case 5
            parameters.colormap='cool';            
            parameters.colormap_label='time';
      end
    
    if get(h.chb_image_load,'Value')==1
        im0=imread(h.fullFileName,1);
    end

    m2=filter_sliders(h,m0);
    set(h.Nparticles_txt,'String',num2str(size(m2,1)));

    x=m2(:,xcol)+get(h.offsetx_slider,'Value');
    y=m2(:,ycol)+get(h.offsety_slider,'Value');

    axes(h.axes1), hold off
    if get(h.chb_image_load,'Value')==1
        im0=imread(h.fullFileName,1);
        if get(h.load_options,'Value')==3
            im0=im0';
        end
%         imagesc(im0),hold,
        imagesc(flipud(rot90(im0,1))),hold,
        set(gca,'YDir','normal')

    end
    if v<4
        
        plot(x,y,'ro','MarkerSize',2)
    else
        figure,
        max_image=1000;
        if get(h.chb_image_load,'Value')==1
            im0=imread(h.fullFileName,1);
            if get(h.load_options,'Value')==3
                im0=im0';
            end
            max_im0=max(max(im0));
            im1=rot90(max_image.*im0./max_im0,-1);
            imagesc(flipud(rot90(im0,1))),hold,
            set(gca,'YDir','normal')

        end
        hold on,
        t=m2(:,1)+max_image*2;
        colormap jet
        scatter(x,y,5,t)
        axis equal
        colorbar
    end
    
        
    xlabel('x,px','color',h.xlabelcolor)
    ylabel('y,px','color',h.xlabelcolor)
    grid
    axis equal
    set(gca,'Xcolor',[0.5 0.5 0.5]);
    set(gca,'Ycolor',[0.5 0.5 0.5]);
else
    status(h,'r','No file loaded!')
end

% Update handles structure
setcallbacks(h);
status(h,'g',strcat('3d data plotted in :',num2str(toc,2),' s'));
