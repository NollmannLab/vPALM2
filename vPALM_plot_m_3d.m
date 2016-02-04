function vPALM_plot_m_3d(h,m0,parameters)
tic

if get(h.chb_loc_by_moment,'Value')==1
    xcol=2;
    ycol=3;
else
    xcol=10;
    ycol=11;
end

if get(h.chb_image_load,'Value')==1
    im0=imread(h.fullFileName,1);
end
    z_code=get(h.zplotoption,'Value');

    if parameters.fit_3dxcorr==1 %      plots 3d xcorr localizations
        m2=filter_sliders(h,m0);
        set(h.Nparticles_txt,'String',num2str(size(m2,1)));
        x=m2(:,xcol);
        y=m2(:,ycol);
        z=m2(:,13);
        zlabelt='Z, nm';
        c=z;
        I=m2(:,4);
        t=m2(:,1);
        res=m2(:,5);
    elseif parameters.fit_3dgaussian==1 %      plots 3d gaussian localizations
        m=filter_3dgauss(m0,parameters,h);
        set(h.Nparticles_txt,'String',num2str(size(m,1)));
        r=m(:,14)./m(:,13);

        if get(h.chb_3dcalapplied,'Value')==0
            m=filter_3dgauss(m0,parameters,h);
            set(h.Nparticles_txt,'String',num2str(size(m,1)));
            zlabelt='Wx/Wy';
            m(:,15)=m(:,14)./m(:,13);
            disp('no calibration file used')
        else % uses cal
            m=filter_3dgauss(h.m_3d,parameters,h);
            zlabelt='Z, nm';
            disp('using z calibration')
        end

        m2=filter_sliders(h,m);
        set(h.Nparticles_txt,'String',num2str(size(m2,1)));
        x=m2(:,xcol);
        y=m2(:,ycol);
        wx=m2(:,13);
        wy=m2(:,14);
        r=wx./wy;
        I=m2(:,4);
        z=m2(:,15);
        t=m2(:,1);
        res=m2(:,5);
    end


    if get(h.chb_plot3d,'Value')==1
        axes(h.axes2), hold off
        switch z_code
            case 1
                c=I;
                parameters.colormap_label='Intensity, a.u.';
                zlabel(parameters.colormap_label,'color',h.xlabelcolor) 

            case 2
                c=z;
                parameters.colormap_label='Z, nm';
                zlabel(parameters.colormap_label,'color',h.xlabelcolor)
            case 3
                density=vPALM_density3d(h,x,y,z);
                parameters.colormap_label='Density, a.u.';
                c=density;
                zlabel(parameters.colormap_label,'color',h.xlabelcolor)
            case 4
                c=res;
                parameters.colormap_label='Resolution, nm';
                zlabel(parameters.colormap_label,'color',h.xlabelcolor) 
            case 5
                c=t;
                parameters.colormap_label='Time, frames';
                zlabel(parameters.colormap_label,'color',h.xlabelcolor) 
        end

        list=get(h.display3D,'String');
        switch list{get(h.display3D,'Value')}
            case 'scatter3'
                scatter3(x,y,z,4,c)
                colormap(parameters.colormap);
                h1 = colorbar;
                ylabel(h1, parameters.colormap_label)
                xlabel('X,px','color',h.xlabelcolor)
                ylabel('Y,px','color',h.xlabelcolor)
                zlabel(zlabelt,'color',h.xlabelcolor)
                set(h.axes2,'Xcolor',[.5 .5 .5]);
                set(h.axes2,'Ycolor',[0.5 0.5 0.5]);
                view(str2num(get(h.AZangle,'String')),str2num(get(h.ELangle,'String')) );
            otherwise
                makeVol_vPALM(x,y,z,str2num(get(h.Nvoxels,'String')),0,1,h);
        end
        
end

% Update handles structure
setcallbacks(h);
status(h,'g',strcat('3d data plotted in :',num2str(toc,2),' s'));
