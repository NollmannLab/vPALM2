function vPALM_stats(hObject, eventdata, h)

list=get(h.stats,'String');
switch list{get(h.stats,'Value')}

    case 'N localizations' %  N local
        status(h,'g','calculating N local statistics');
        h.m2=filter_sliders(h,h.m);

        axes(h.axes2)
        hold off,
        [n,x]=hist(h.m2(:,1),size(h.m2(:,1),1));
        bar(x, n)
        xlabel('Frame #','color',h.xlabelcolor)
        ylabel('PDF (a.u.)','color',h.xlabelcolor)
        grid
        set(gca,'Xcolor',[0.5 0.5 0.5]);
        set(gca,'Ycolor',[0.5 0.5 0.5]);


    case 'N local/pixel'% N local/ pixel
        status(h,'g','calculating N local/px statistics...');
        h.m2=filter_sliders(h,h.m);

        axes(h.axes2)
        hold off,
        m=filter_sliders(h,h.m2);
        x=m(:,2); y=m(:,3);
%         det_thresh=580/2/1.4./str2num(get(h.pixelsize,'String'));
        det_thresh=str2num(get(h.h_sliders.detect_thresh,'String'))./str2num(get(h.pixelsize,'String'));

%         det_thresh=  str2num(get(h.detect_thresh,'String'))./str2num(get(h.pixelsize,'String'));
%         tic
        for iparticle=1:size(x)
            nfr= find ( m(:,1)== m(iparticle,1) ); % particles in same frame
            n= find( sqrt( (x(iparticle)-x(nfr)).^2 + (y(iparticle)-y(nfr)).^2 ) < 2*det_thresh);  %earticles in same PSF    
            density(iparticle)=size(n,1 );
        end
%         toc
        
        [n,x]=hist(density);
        bar(x, n);
        xlabel('Frame #','color',h.xlabelcolor)
        ylabel('PDF (a.u.)','color',h.xlabelcolor)
        grid
        set(gca,'Xcolor',[0.5 0.5 0.5]);
        set(gca,'Ycolor',[0.5 0.5 0.5]);
        disp('done');
    case 'density' %density
        status(h,'g','calculating density statistics');
        h.m2=filter_sliders(h,h.m);

        axes(h.axes2)
        hold off,
        m=filter_sliders(h,h.m2);
        x=m(:,2); y=m(:,3);
        det_thresh=580/2/1.4./str2num(get(h.pixelsize,'String'));
%         str2num(get(h.detect_thresh,'String'))./str2num(get(h.pixelsize,'String'));

        % N particles in a sphere with a radius defined by det_thresh (PSF)
        
        for iparticle=1:size(x)
            n= find( sqrt( (x(iparticle)-x).^2 + (y(iparticle)-y).^2 ) < det_thresh);
            
            density(iparticle)=size(n,1 );
        end

        [n,x]=hist(density);
        bar(x, n)
        xlabel('# events within a PSF','color',h.xlabelcolor)
        ylabel('PDF (a.u.)','color',h.xlabelcolor)
        grid
        set(gca,'Xcolor',[0.5 0.5 0.5]);
        set(gca,'Ycolor',[0.5 0.5 0.5]);
    
    case 'dy' %dy
        status(h,'g','calculating dy statistics');
        h.m2=filter_sliders(h,h.m);axes(h.axes2)
        hold off,
        m=filter_sliders(h,h.m2);
        y=m(:,3);
        hist((y-mean(y))*100)
        disp('dy');
        xlabel('y frequency, px','color',h.xlabelcolor)
        ylabel('PDF (a.u.)','color',h.xlabelcolor)
        grid
        set(gca,'Xcolor',[0.5 0.5 0.5]);
        set(gca,'Ycolor',[0.5 0.5 0.5]);
        
        axes(h.axes1)
        ynm=y*str2num(get(h.pixelsize,'String'));
        hold off,
        plot(y)
        hold on,
        plot(meanfilter(100,y),'k','Linewidth',2)
        xlabel('frame #','color',h.xlabelcolor)
        ylabel('y, px','color',h.xlabelcolor)
        grid
        set(gca,'Xcolor',[0.5 0.5 0.5]);
        set(gca,'Ycolor',[0.5 0.5 0.5]);
        disp(strcat('drift corrected FWHM_x=',num2str(std(ynm-meanfilter(100,ynm)'))));
        disp(strcat('average drift x (nm)=',num2str(2.2*std(ynm))));


        status(h,'g',strcat('FWHM_y=',num2str(2.2*100*std(y-mean(y))),' nm'));
    case 'dx' % dx
        status(h,'g','calculating dx statistics');
        h.m2=filter_sliders(h,h.m);axes(h.axes2)
        hold off,
        m=filter_sliders(h,h.m2);
        x=m(:,2); 
        hist((x-mean(x))*100)
        disp('dx');
         xlabel('x frequency, px','color',h.xlabelcolor)
        ylabel('PDF (a.u.)','color',h.xlabelcolor)
        grid
        set(gca,'Xcolor',[0.5 0.5 0.5]);
        set(gca,'Ycolor',[0.5 0.5 0.5]);
        
        axes(h.axes1)
        xnm=x*str2num(get(h.pixelsize,'String'));
        hold off,
        plot(x)
        hold on,
        plot(meanfilter(100,x),'k','Linewidth',2)
        xlabel('frame #','color',h.xlabelcolor)
        ylabel('x, px','color',h.xlabelcolor)
        grid
        set(gca,'Xcolor',[0.5 0.5 0.5]);
        set(gca,'Ycolor',[0.5 0.5 0.5]);
        disp(strcat('drift corrected FWHM_x=',num2str(std(xnm-meanfilter(100,xnm)'))));
        disp(strcat('average drift x (nm)=',num2str(2.2*std(xnm))));

        status(h,'g',strcat('FWHM_x=',num2str(2.2*100*std(x-mean(x))),' nm'));
        
    case 'dz' % dx
        status(h,'g','calculating dz statistics');
        h.m2=filter_sliders(h,h.m);
        axes(h.axes2)
        hold off,
        m=filter_sliders(h,h.m2);
        z=m(:,15); 
        hist((z-mean(z)))
        disp('dz');
        xlabel('z frequency, nm','color',h.xlabelcolor)
        ylabel('PDF (a.u.)','color',h.xlabelcolor)
        grid
        set(gca,'Xcolor',[0.5 0.5 0.5]);
        set(gca,'Ycolor',[0.5 0.5 0.5]);
        
        axes(h.axes1)
        hold off,
        plot(z)
        hold on,
        plot(meanfilter(100,z),'k','Linewidth',2)
%         plot(z-meanfilter(100,z)','r','Linewidth',2)
        xlabel('frame #','color',h.xlabelcolor)
        ylabel('z, nm','color',h.xlabelcolor)
        grid
        set(gca,'Xcolor',[0.5 0.5 0.5]);
        set(gca,'Ycolor',[0.5 0.5 0.5]);
                
                
        status(h,'g',strcat('FWHM_z=',num2str(std(z)),' nm'));
        frsz=[];
        for i=1:100:size(z,1)
            if i+100-1>size(z,1)
                range=[i:1:size(z,1)];
            else 
                range=linspace(i,i+100-1,100);
            end
            frsz=[frsz,std(z(range))];
        end
        
        disp(strcat('fractional FWHM_z=',num2str(mean(frsz))));
        disp(strcat('drift corrected FWHM_z=',num2str(std(z-meanfilter(100,z)'))));

                
    case 'resolution' %'resolution'
        status(h,'g','calculating dy statistics');
        h.m2=filter_sliders(h,h.m);axes(h.axes2)
        hold off,
        m=filter_sliders(h,h.m2);
        res=m(:,5);
        res=res(find(res<250));
        [hx,b]=hist(res,20);
        bar(b,hx)
        disp('resolution');
        xlabel('resolution, nm','color',h.xlabelcolor)
        ylabel('PDF (a.u.)','color',h.xlabelcolor)
        grid
        axis([0 200 0 max(hx)*1.2])

        set(gca,'Xcolor',[0.5 0.5 0.5]);
        set(gca,'Ycolor',[0.5 0.5 0.5]);
%         status(h,'g',strcat('mean res',num2str(mean(res)),' nm'));        
    case 'background' %,'background'
        status(h,'g','calculating dy statistics');
        h.m2=filter_sliders(h,h.m);axes(h.axes2)
        hold off,
        m=filter_sliders(h,h.m2);
        res=m(:,6);
        hist((res))
        disp('background');
        xlabel('background','color',h.xlabelcolor)
        ylabel('PDF (a.u.)','color',h.xlabelcolor)
        grid
        set(gca,'Xcolor',[0.5 0.5 0.5]);
        set(gca,'Ycolor',[0.5 0.5 0.5]);
        status(h,'g',strcat('mean res',num2str(mean(res)),' nm'));        
               
        
    case 'Intensity' % ,'Intensity'
        status(h,'g','calculating dy statistics');
        h.m2=filter_sliders(h,h.m);axes(h.axes2)
        hold off,
        m=filter_sliders(h,h.m2);
        maxt=max(h.m(:,1));
        loc_times=m(:,1);
        t=[1:maxt];
        intensity=zeros(1,size(t,2));

        for i=1:size(loc_times,1)
            intensity(loc_times(i))=m(i,4);
        end

        plot(t,intensity)
        disp('Intensity');
        xlabel('Intensity','color',h.xlabelcolor)
        ylabel('PDF (a.u.)','color',h.xlabelcolor)
        grid
        axis([0 maxt 0 max(intensity)]);
        set(gca,'Xcolor',[0.5 0.5 0.5]);
        set(gca,'Ycolor',[0.5 0.5 0.5]);
        status(h,'g',strcat('mean intensity',num2str(mean(intensity)),' a.u.'));        
               
        
end

