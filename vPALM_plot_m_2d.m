function vPALM_plot_m_2d(h,m,parameters)
tic
    z_code=get(h.zplotoption,'Value');
    
    status(h,'r','working...')
    m2=filter_sliders(h,m);
    set(h.Nparticles_txt,'String',num2str(size(m2,1)));

    x=m2(:,2)+get(h.offsetx_slider,'Value');
    y=m2(:,3)+get(h.offsety_slider,'Value');
    
%     axes(h.axes1), hold off
%     if get(h.chb_image_load,'Value')==1
%         im0=imread(h.fullFileName,1);
%         if get(h.load_options,'Value')==3
%             im0=im0';
%         end
%         imagesc(im0),hold,
%     end
%     
%     plot(x,y,'r.','MarkerSize',4)
%     xlabel('x,px','color',[.8 .8 .8])
%     ylabel('y,px','color',[.8 .8 .8])
%     grid
%     set(gca,'Xcolor',[0.5 0.5 0.5]);
%     set(gca,'Ycolor',[0.5 0.5 0.5]);

    if get(h.chb_plot3d,'Value')==1
        axes(h.axes2), hold off

        switch z_code
            case 1
                z=m2(:,4);
                c=z;
                zlabel('Intensity, a.u.','color',h.xlabelcolor) ;
                parameters.colormap_label='Intensity, a.u.';

            case 2
                z=ones(size(m2(:,4)));
                c=z;
                parameters.colormap_label='N/A';
                zlabel('none','color',h.xlabelcolor) ;
            case 3
                density=density2d(h,x,y);
                z=density;
                c=z;
                parameters.colormap_label='Density, a.u.';
                zlabel('Density, a.u.','color',h.xlabelcolor) ;
            case 4
                z=m2(:,5);
                parameters.colormap_label='Resolution, nm';
                c=z;
                zlabel('Resolution, nm','color',h.xlabelcolor)
            case 5
                z=m2(:,1);
                c=z;
                parameters.colormap_label='Time, frames';
                zlabel('time, frames','color',h.xlabelcolor)
        end
        

        if get(h.chb_voxels,'Value')==1
            makeVol_vPALM(x,y,z,str2num(get(h.Nvoxels,'String')),0,0,h);
        else
            scatter3(x,y,z,4,c)
            xlabel('X,px','color',h.xlabelcolor)
            ylabel('Y,px','color',h.xlabelcolor)
            zlabel(parameters.colormap_label,'color',h.xlabelcolor)
            colormap(parameters.colormap);
            h1 = colorbar;
            ylabel(h1, parameters.colormap_label)
            set(h.axes2,'Xcolor',[0.5 0.5 0.5]);
            set(h.axes2,'Ycolor',[0.5 0.5 0.5]);
            view(str2num(get(h.AZangle,'String')),str2num(get(h.ELangle,'String')) );
        end
    end
    
    status(h,'g',strcat('2d data plotted in :',num2str(toc,2),' s'));

%     colormap parameters.colormap;
