% vPALM GUI
function varargout = vPALM2(varargin)
% vPALM Brief description of GUI.
%       Comments displayed at the command line in response 
%       to the help command. 

% (Leave a blank line following the help.)
% yet to do:

% slider to change
% add 2d plots with particles detected alone on right panel
% option menu to select color code in 2d or 3d plots
%   number of molecules within given radius (clustering)
%   z
%   Intensity of localization
%   eccentricity in 3d
%   
% variable to select the 'clustering radius' (sort of density plot)
%
% 02012013- to do:
% 1- popup with options for loading from Mtt, qpalm, and ens.
% 2- ROI in axes1 to select ROI in xy
% 4- additional pushbutton to open image from file independently.
% 4- modify qPALM so that final values after calibration are passed in the
% m matrix (moslty when using the 3d gaussian fitting)
% 5- sliders to select density (think about it first)

%%  Initialization tasks

% mPropertyDefs = {...        
%                 'iconwidth',  @localValidateInput, 'mIconWidth';
%                 'iconheight', @localValidateInput, 'mIconHeight';
%                 'iconfile',   @localValidateInput, 'mIconFile'};
% mIconWidth = 16;   % Use input property 'iconwidth' to initialize
% mIconHeight = 16;  % Use input property 'iconheight' to initialize
% mIconFile = fullfile(matlabroot,'toolbox/matlab/icons/'); 
                   % Use input property 'iconfile' to initialize
scrsz = get(0,'ScreenSize');
% FIGH=scrsz(4)/1;
% FIGW=scrsz(3)*.8;
FIGH=746;
FIGW=800; %1092


h.MainFig = figure(...
              'Name','vPALM',...
              'Units','pixels',...
              'MenuBar','none',...
              'Toolbar','figure',...
              'Color','w',...              
              'Position',[1 1 FIGW FIGH],...
              'Visible','off');
          

% f = figure('Visible','off','Name','My GUI',...
%            'Position',[360,500,450,285]);

%%  Construct the components
h=vPALM2_gui_construction(h,FIGH,FIGW);

%%  Initialization tasks

% setcallbacks(h);

set(h.MainFig,'Visible','on')

set(h.MainFig, 'UserData', h);

%  Callbacks for MYGUI

h.h_vis=vPALM2_2d3d_visualization(h);

h.h_sliders=vPALM2_sliders(h);

set(h.MainFig,'CloseRequestFcn',{@vPALM2_newclosefunction,h});

% set(h.h_sliders.detect_thresh,'String','-10');

%% ends
setcallbacks(h);

disp('GUI initialisation OK');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vPALM2_newclosefunction(hObject, eventdata,h)
disp('everything is closing down!');

delete(h.h_vis.fig);
delete(h.h_sliders.fig);
delete(h.MainFigure3D);
delete(h.MainFig);

%% load_file
function h=load_file(hObject,eventdata,h)

disp_nice('Give the filename of the mat file')
% write_parameters(hObject,h,parameters);

MAX_r0=10;
MAX_ELLIP=5;

set(h.chb_3dcalapplied,'Value',0);
set(h.chb_drift_correction,'Value',0);

list=get(h.load_options,'String');
switch list{get(h.load_options,'Value')}
    case 'vPALM' %vPALM
        
        [fullFileName,path1] = uigetfile('*vPALM.mat','Open MAT File');
        cd(path1);
        load(fullFileName);
        set(h.chb_drift_correction,'Value',parameters.chb_drift_correction);
        set(h.chb_image_load,'Value',parameters.chb_image_load);
        set(h.pixelsize,'String',num2str(parameters.pixelsize));
        set(h.chb_cal_3d,'Value',parameters.chb_cal_3d);
        set(h.chb_3dcalapplied,'Value',parameters.chb_3dcalapplied);

        clear_m_matrices;

        h.m=[];
  
        h.m(:,1)=data.Nframe; % frame number
        h.m(:,2)=data.x./str2num(get(h.pixelsize,'String'))+1; % pos x
        h.m(:,3)=data.y./str2num(get(h.pixelsize,'String'))+1; % pos y
        h.m(:,4)=data.Intensity; % intensity
        h.m(:,5)=data.Resolution; % resolution, nm
        h.m(:,6)=data.Backround; %background
        h.m(:,7)=data.Molecule_id; %molecule
        h.m(:,8)=data.Molecule2; %molecule
        h.m(:,9)=data.Channel; %channel
        h.m(:,10)=data.Slice; %slice
        h.m(:,12)=data.A; %a ???
        
        if parameters.chb_3dcalapplied==0
            h.m(:,11)=data.Width; %width
            h.m(:,13)=data.Theta; %theta
            h.m(:,14)=data.Xpos; %xpos
            h.m(:,15)=data.Ypos; %ypos
        else
            h.m(:,11)=data.Asymmetry; %asymetry
            h.m(:,13)=data.Wx; %wx
            h.m(:,14)=data.Wy; %wy
            h.m(:,15)=data.z; %z    
            h.parameters.v_gauss_cal=1;
            h.parameters.min_ellip=0;
            h.parameters.max_ellip=MAX_ELLIP;
            h.parameters.min_w=0 ;
            h.parameters.max_w=MAX_r0*str2num(get(h.pixelsize,'String'));


        end
        
        h.fullFileName=parameters.fullFileName;
        h.parameters.fit_2dgaussian=parameters.fit_2dgaussian;
        h.parameters.fit_3dgaussian=parameters.fit_3dgaussian;
        h.parameters.fit_3dxcorr=parameters.fit_3dxcorr;
        h.parameters.fullFileName=parameters.fullFileName;
        h.parameters.pixelsize=parameters.pixelsize;
        h.parameters.chb_drift_correction=parameters.chb_drift_correction;
        h.parameters.chb_image_load=parameters.chb_image_load;
        h.parameters.chb_cal_3d=parameters.chb_cal_3d;
        h.parameters.chb_3dcalapplied=parameters.chb_3dcalapplied;
        h.parameters.v_gauss_cal=parameters.v_gauss_cal;


        if parameters.chb_3dcalapplied==1 % loads calibration parameters if 3dcal was available upon saving
            h.wxPoly=parameters.wxPoly ;
            h.wyPoly=parameters.wyPoly  ;
            h.ccPolyWtoZ=parameters.ccPolyWtoZ;
            h.m_3d=h.m;
        end
        
        if get(h.chb_drift_correction,'Value')==0
            ;
        else % both m-matrices will be the same to avoid crashing if one is not there.
            h.m_drift=h.m;
        end
        
        
        

    case 'qPALM' %qPALM
        [fullFileName,path1] = uigetfile('*_analysed.mat','Open MAT File');
        cd(path1);
        load(fullFileName);
        h.m=matrice_results';
        h.fullFileName=parameters.fullFileName;
        h.parameters=parameters;
        
        if parameters.v_gauss_cal==1
          set(h.chb_cal_3d,'Value',1);
        else
          set(h.chb_cal_3d,'Value',0);
        end
    case 'PALMcbs' % PALMcbs
        [fullFileName,path1] = uigetfile('*_analysed.mat','Open MAT File');
        cd(path1);
        load(fullFileName);
        h.m=AllData.m1';
        h.m(:,2)=AllData.m1(3,:);
        h.m(:,3)=AllData.m1(2,:);

        h.fullFileName=fullFileName;
        parameters.fit_2dgaussian=1;
        parameters.fit_3dgaussian=0;
        parameters.fit_3dxcorr=0;
        h.parameters=parameters;
        if get(h.chb_image_load,'Value')==1
            disp_nice('Give the filename of the image file')
            [TIFfilename,path1] = uigetfile('*.tif','Open TIF File');
            h.fullFileName=TIFfilename;
            parameters.fullFileName=TIFfilename;
        else
            h.fullFileName=fullFileName; 
            parameters.fullFileName=fullFileName;
        end

    
    case 'MTT' %4 %MTT
        [fullFileName,path1] = uigetfile('*_analysed.mat','Open MAT File');
        cd(path1);
        load(fullFileName);
        clear_m_matrices;

        switch size(matrice_results,1)
            case 5
                % frame, x, y, alpha, r0
                parameters.fit_2dgaussian=1;
                parameters.fit_3dgaussian=0;
                parameters.fit_3dxcorr=0; 
                m=matrice_results';

            case 8
                % frame, x, y, alpha, 0, r0, rx, ry 
                
                mp=matrice_results';
                clear m;
                correct_range=find(mp(:,6)<MAX_r0); % filters particles with r0>5 pixels
                
                m0=mp(correct_range,:);
                m=m0;
                m(:,5)=m0(:,6); % r0, nm 
                m(:,6)=m0(:,5); %background
                m(:,13)=m0(:,7); %Wx %%added 
                m(:,14)=m0(:,8); %Wy %%added 
                m(:,10)=m0(:,5); %angle
                m(:,11)=m0(:,8)./m0(:,7); %asymmetry
                m(:,12)=m0(:,5); %channel
                m(:,8)=m0(:,5); %slice
                m(:,9)=m0(:,5); %pos

                frames=m(:,1);
                [~,index_m]=sort(frames);
                m2(:,:)=m(index_m,:);

                h.m=m2;
                parameters.fit_2dgaussian=0;
                parameters.fit_3dgaussian=1;
                parameters.fit_3dxcorr=0;
                
                parameters.v_gauss_cal=1;
                parameters.min_ellip=0;
                parameters.max_ellip=MAX_ELLIP;
                parameters.min_w=0 ;
                parameters.max_w=MAX_r0*str2num(get(h.pixelsize,'String'));

        end
            
        if get(h.chb_collate,'Value')==0
            h.m=m; % does not collate
        else
            % collates files
            MaxFrame = max(h.m(:,1));
            m(:,1) = m(:,1) + MaxFrame;
            h.m=[h.m;m];
        end
        
        if get(h.chb_image_load,'Value')==1
            disp_nice('Give the filename of the image file')
            [TIFfilename,path1] = uigetfile('*.tif','Open TIF File');
            h.fullFileName=TIFfilename;
            parameters.fullFileName=TIFfilename;
        else
            h.fullFileName=fullFileName; 
            parameters.fullFileName=fullFileName;
        end
        
        
        h.parameters=parameters;


    case 'rapidStorm' %rapidSTORM
        [fullFileName,path1] = uigetfile('*.*','Open rapidstorm File');
        cd(path1);
        A=importdata(fullFileName);
        
        clear_m_matrices;

        %1,2 x and y in nm
        %3 position in sample space in z
        %4 frame number
        %5 emission strength" unit="A/D count
        %6 PSF FWHM in x dimension" unit="nanometer"
        %7 PSF FWHM in x dimension" unit="nanometer"
        %8 fit residue chi square value" unit="dimensionless
        %9 local background" unit="A/D count
  
        m(:,1)=A.data(:,4); % frame number
        m(:,2)=1+A.data(:,1)./str2num(get(h.pixelsize,'String')); % pos x in px
        m(:,3)=1+A.data(:,2)./str2num(get(h.pixelsize,'String')); %pos y in px
        m(:,4)=A.data(:,5); % intensity
        m(:,5)=zeros(size(A.data(:,1))); % resolution, nm
        m(:,6)=A.data(:,9); %background
        m(:,13)=A.data(:,6)./str2num(get(h.pixelsize,'String')); %Wx %%added tue1405
        m(:,14)=A.data(:,7)./str2num(get(h.pixelsize,'String')); %Wy %%added tue1405
        m(:,10)=A.data(:,4); %angle
        m(:,11)=A.data(:,6)./A.data(:,7); %asymmetry
        m(:,12)=ones(size(A.data(:,1))); %channel
        m(:,8)=ones(size(A.data(:,1))); %slice
        m(:,9)=zeros(size(A.data(:,1))); %pos
        m(:,7)=A.data(:,8); % fit residue, chi square
        
        frames=m(:,1);
        [~,index_m]=sort(frames);
        m2(:,:)=m(index_m,:);
        
        if get(h.chb_collate,'Value')==0
            h.m=m2; % does not collate
        else
            % collates files
            MaxFrame = max(h.m(:,1));
            m2(:,1) = m2(:,1) + MaxFrame;
            h.m=[h.m;m2];
        end
        
        h.fullFileName=fullFileName;
        h.parameters.fit_2dgaussian=0;
        h.parameters.fit_3dgaussian=1;
        h.parameters.fit_3dxcorr=0;
        h.parameters.v_gauss_cal=1;
        h.parameters.min_ellip=0;
        h.parameters.max_ellip=MAX_ELLIP;
        h.parameters.min_w=0 ;
        h.parameters.max_w=MAX_r0*str2num(get(h.pixelsize,'String'));

        h.parameters.fullFileName=fullFileName;
        set(h.chb_image_load,'Value',0);
        
        
    case 'micromanager-2D'%6 %micromanager 2d
        
        [fullFileName,path1] = uigetfile('*.*','Open DAT File');
        cd(path1);
        A=importdata(fullFileName);
        clear m; 

        m(:,1)=A.data(:,4); % frame number
        m(:,2)=1+A.data(:,6)./str2num(get(h.pixelsize,'String')); % pos x
        m(:,3)=1+A.data(:,7)./str2num(get(h.pixelsize,'String')); %pos y
        m(:,4)=A.data(:,8); % intensity
        m(:,5)=A.data(:,15); % resolution, nm
        m(:,6)=A.data(:,9); %background
        m(:,7)=A.data(:,1); %molecule
        m(:,8)=A.data(:,2); %molecule
        m(:,9)=A.data(:,3); %channel
        m(:,10)=A.data(:,5); %slice
        m(:,11)=A.data(:,10); %width
        m(:,12)=A.data(:,11); %a
        m(:,13)=A.data(:,12); %theta
        m(:,14)=A.data(:,13); %xpos
        m(:,15)=A.data(:,14); %ypos

%         m(:,7:14)=[zeros(size(m,1),6),m(:,4),zeros(size(m,1),1)];
        
        if get(h.chb_collate,'Value')==0
            h.m=m;% does not collate
        else
            % collates files
            MaxFrame = max(h.m(:,1));
            m(:,1) = m(:,1) + MaxFrame;
            h.m=[h.m;m];
        end
        
        h.fullFileName=fullFileName;
        h.parameters.fit_2dgaussian=1;
        h.parameters.fit_3dgaussian=0;
        h.parameters.fit_3dxcorr=0;
        h.parameters.fullFileName=fullFileName;
%         h.parameters=parameters;
        set(h.chb_image_load,'Value',0);
        
 case 'micromanager-3D'% 7 %micromanager 3d
        [fullFileName,path1] = uigetfile('*.*','Open DAT File');
        cd(path1);
        A=importdata(fullFileName)
        clear m;
        m(:,1)=A.data(:,4); % frame number
        m(:,2)=1+A.data(:,6)./str2num(get(h.pixelsize,'String')); % pos x in px
        m(:,3)=1+A.data(:,7)./str2num(get(h.pixelsize,'String')); %pos y in px
        m(:,4)=A.data(:,8); % intensity
        m(:,5)=A.data(:,15); % resolution, nm
        m(:,6)=A.data(:,9); %background
        m(:,13)=A.data(:,10)./str2num(get(h.pixelsize,'String')); %Wx %%added tue1405
        m(:,14)=A.data(:,10)./A.data(:,11)./str2num(get(h.pixelsize,'String')); %Wy %%added tue1405
        m(:,10)=A.data(:,12); %angle
        m(:,11)=A.data(:,11); %asymmetry
        m(:,12)=A.data(:,2); %channel
        m(:,8)=A.data(:,3); %slice
        m(:,9)=A.data(:,5); %pos
        
        frames=m(:,1);
        [~,index_m]=sort(frames);
        m2(:,:)=m(index_m,:);
        
        
        if get(h.chb_collate,'Value')==0
            h.m=m2;
        else
            % collates files
            MaxFrame = max(h.m(:,1));
            m2(:,1) = m2(:,1) + MaxFrame;
            h.m=[h.m;m2];
        end
        
        h.fullFileName=fullFileName;
        h.parameters.fit_2dgaussian=0;
        h.parameters.fit_3dgaussian=1;
        h.parameters.fit_3dxcorr=0;
        h.parameters.v_gauss_cal=1;
        h.parameters.min_ellip=0;
        h.parameters.max_ellip=MAX_ELLIP;
        h.parameters.min_w=0 ;
        h.parameters.max_w=MAX_r0*str2num(get(h.pixelsize,'String'));

        h.parameters.fullFileName=fullFileName;
%         h.parameters=parameters;
        set(h.chb_image_load,'Value',0);

end

        
status(h,'g','everything went apparently fine')

set(h.filename,'String',fullFileName,'ForegroundColor',[.7 .95 .7]);
set(h.chb_drift_correction,'Value',0);

set(h.Nparticles_txt,'String',num2str(size(h.m,1)));

set(h.chb_data_loaded,'Value',1);

slider_setup(h,h.m, h.parameters);

       
% Update handles structure
setcallbacks(h)


%% clear m_matrices
function clear_m_matrices()

        if exist('h.m')==1
            clear h.m;
        end
        
        if exist('h.m_drift')==1
            clear h.m_drift;
        end
        
        if exist('h.m2')==1
            clear h.m2;
        end

        if exist('m')==1
            clear m;
        end
        

%% plot_localizations
function h=plot_localizations_3d(hObject,eventdata,h)
    status(h,'r','working...')

if get(h.chb_data_loaded,'Value')==1
    
    m0=h.m;
    parameters=h.parameters;
    v=get(h.Zcolormap,'Value');
    switch v
        case 0
            parameters.colormap='jet';
        case 1
            parameters.colormap='hot';
        case 2
            parameters.colormap='hsv';
        case 3
            parameters.colormap='gray';
    end
    
    if parameters.fit_2dgaussian==1 | (parameters.fit_3dgaussian==0 & parameters.fit_3dxcorr==0)
        plot_m_2d(h,m0,parameters)
    elseif (parameters.fit_3dgaussian==1 | parameters.fit_3dxcorr==1)
        plot_m_3d(h,m0,parameters);
    end

else
    status(h,'r','No file loaded!')
end

% Update handles structure
setcallbacks(h)

%% plot_localizations_2d
function h=plot_localizations_2d(hObject,eventdata,h)
tic
status(h,'r','working...');

if get(h.chb_data_loaded,'Value')==1
    
    m0=h.m;
    parameters=h.parameters;
    v=get(h.Zcolormap,'Value');
    switch v
        case 0
            parameters.colormap='jet';
        case 1
            parameters.colormap='hot';
        case 2
            parameters.colormap='hsv';
        case 3
            parameters.colormap='gray';
    end
    
    if get(h.chb_image_load,'Value')==1
        im0=imread(h.fullFileName,1);
    end

    m2=filter_sliders(h,m0);
    set(h.Nparticles_txt,'String',num2str(size(m2,1)));

    x=m2(:,2)+get(h.offsetx_slider,'Value');
    y=m2(:,3)+get(h.offsety_slider,'Value');
    
    axes(h.axes1), hold off
    if get(h.chb_image_load,'Value')==1
        im0=imread(h.fullFileName,1);
        if get(h.load_options,'Value')==3
            im0=im0';
        end
        imagesc(im0),hold,
    end
    plot(x,y,'ro','MarkerSize',2)
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

%%
function cal_3d(hObject, eventdata, h)

tic
[h.wxPoly,h.wyPoly,h.ccPolyWtoZ]=vPALM2_cal3D(h);

%     m=h.m;
%     set(h.Nparticles_txt,'String',num2str(size(m,1)));
% 
%     r=m(:,11);
%     wx=m(:,13);
%     wy=m(:,14);
%     fr=m(:,1);
%     
%     figure, plot(fr,wx,'or-')  , hold on, plot(fr,wy,'ob-')    
%     disp('Click first point of event to analyze:')
%     regionstart= ginput(1);
%     regionstart= find(fr> regionstart(1) );
%     regionstart= min(regionstart);
%     regionstart_Z= fr(regionstart);
%     disp('Click last point of event to analyze:')
%     regionend= ginput(1);
%     regionend= find(fr < regionend(1) );
%     regionend= max(regionend);
%     regionend_Z= fr(regionend);   
%     
%     range=regionstart:regionend;
%     stepsize=str2num(get(h.stepsize,'String'));
%     zcal=fr(range);
%     xcal=wx(range);
%     ycal=wy(range);
%     
%     [ax,p]=polyfit(zcal,xcal,5);
%     plot(zcal,polyval(ax,zcal),'b','Linewidth',5);
%     
%     [ay,p]=polyfit(zcal,ycal,5);
%     plot(zcal,polyval(ay,zcal),'r','Linewidth',5);
%     
%     xlabel('piezo step')
%     ylabel('Wx, Wy, px')

%% display
%         zgrid=stepsize*linspace(min(zcal),max(zcal),1000);
%         wx_cal=polyval(ax,zgrid);
%         wy_cal=polyval(ay,zgrid);
%         minsq2=[];
%         for ievent=1:size(m,1)
%             [~,minsq2(ievent)]=min((wx(ievent)-wx_cal).^2+(wy(ievent)-wy_cal).^2);
%         end
%         z=minsq2*(parameters.cal_to-parameters.cal_from)/1000*parameters.stepsize;
%     h.zcal_xpars=ax;
%     h.zcal_ypars=ay; 
%     h.parameters.v_gauss_cal=1;
%     h.parameters.stepsize=str2num(get(h.stepsize,'String'));
% %     h.parameters.cal_from=0;
% %     h.parameters.cal_to=size(zcal,1)*str2num(get(h.stepsize,'String'));
%     h.parameters.cal_from=regionstart_Z;
%     h.parameters.cal_to=regionend_Z;
%   
%     z=calc_z_from_cal(h,m,h.parameters);
%     figure, plot(z), 
%     xlabel('rel fr number')
%     ylabel('z position, nm')


    set(h.chb_cal_3d,'Value',1);
    
setcallbacks(h);
status(h,'g',strcat('3d data plotted in :',num2str(toc,2),' s'));

   disp('') ;
    
function plot_m_3d(h,m0,parameters)
tic
    if get(h.chb_image_load,'Value')==1
        im0=imread(h.fullFileName,1);
    end
    z_code=get(h.zplotoption,'Value');

if parameters.fit_3dxcorr==1 %      plots 3d xcorr localizations
    m2=filter_sliders(h,m0);
    set(h.Nparticles_txt,'String',num2str(size(m2,1)));

%     x=m2(:,2);
%     y=m2(:,3);
    x=m2(:,2);
    y=m2(:,3);
    z=m2(:,13);
    zlabelt='Z, nm';
    c=z;
    I=m2(:,4);
    t=m2(:,1);

elseif parameters.fit_3dgaussian==1 %      plots 3d gaussian localizations
    
    m=filter_3dgauss(m0,parameters,h);
    set(h.Nparticles_txt,'String',num2str(size(m,1)));

    r=m(:,14)./m(:,13);

    if get(h.chb_3dcalapplied,'Value')==0
    % if parameters.v_gauss_cal==0 % no calibation file, wil show ratio wy/wx;
%         c=r;
        m=filter_3dgauss(m0,parameters,h);
        set(h.Nparticles_txt,'String',num2str(size(m,1)));

        zlabelt='Wx/Wy';
        m(:,15)=m(:,14)./m(:,13);
%         z=r;
        disp('no calibration file used')
    else % uses cal
%         z=calc_z_from_cal(h,m,h.parameters);
        m=filter_3dgauss(h.m_3d,parameters,h);

        zlabelt='Z, nm';
%         m(:,15)=z;
        disp('using z calibration')
    end
    
    m2=filter_sliders(h,m);
    set(h.Nparticles_txt,'String',num2str(size(m2,1)));

%     if get(h.load_options,'Value')>4
        x=m2(:,2);
        y=m2(:,3);
%     else
%         x=m2(:,10);
%         y=m2(:,11);
%     end
%     
    wx=m2(:,13);
    wy=m2(:,14);
    r=wx./wy;
    I=m2(:,4);
    z=m2(:,15);
    t=m2(:,1);
    
end

    
%     axes(h.axes1), hold off
%     if get(h.chb_image_load,'Value')==1
%         imagesc(im0);
%         hold on
%     else
%         hold off
%     end
% %     scatter3(x,y,z,4,z)
%     plot(x,y,'ro','MarkerSize',2)
% %     colormap parameters.colormap;
%     set(gca,'Xcolor',[0.5 0.5 0.5]);
%     set(gca,'Ycolor',[0.5 0.5 0.5]);
%     xlabel('X,px','color',[.8 .8 .8])
%     ylabel('Y,px','color',[.8 .8 .8])
    
    if get(h.chb_plot3d,'Value')==1

        axes(h.axes2), hold off
        switch z_code
            case 1
                c=I;
                zlabel('Intensity, a.u.','color',h.xlabelcolor) 

            case 2
                c=z;
                zlabel('z, nm','color',h.xlabelcolor)
            case 3
                density=density3d(h,x,y,z);
                c=density;
                zlabel('Density, a.u.','color',h.xlabelcolor)
            case 4
                c=r;
                zlabel('Ellipticity','color',h.xlabelcolor) 
            case 5
                c=t;
                zlabel('Time, fr','color',h.xlabelcolor) 
        end
        
        
        list=get(h.display3D,'String');
        switch list{get(h.display3D,'Value')}
            case 'scatter3'
                scatter3(x,y,z,4,c)
                colorbar
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

%%    
function z=calc_z_from_cal(h,m,parameters);


        ax=h.zcal_xpars;
        ay=h.zcal_ypars;
        zgrid=linspace(parameters.cal_from,parameters.cal_to,1000);
        wx_cal=polyval(ax,zgrid);
        wy_cal=polyval(ay,zgrid);
        wx=m(:,13);
        wy=m(:,14);
      
        
        minsq2=[];
        for ievent=1:size(m,1)
            [min2(ievent),minsq2(ievent)]=min( (wx(ievent)-wx_cal).^2+(wy(ievent)-wy_cal).^2 );
        end
        z=(zgrid(minsq2)-parameters.cal_from)*str2num(get(h.stepsize,'String'));
        z=z-min(z);
%         c=z;
disp('done')

function plot_m_2d(h,m,parametersd)
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

            case 2
                z=ones(size(m2(:,4)));
                c=z;
                zlabel('none','color',h.xlabelcolor) ;
            case 3
                density=density2d(h,x,y);
                z=density;
                c=z;
                zlabel('Intensity, a.u.','color',h.xlabelcolor) ;
            case 4
               z=m2(:,1);
                density=density2d(h,x,y);
                c=density;
                zlabel('time, frames','color',h.xlabelcolor)
            case 5
                z=m2(:,1);
                density=density2d(h,x,y);
                c=z;
                zlabel('time, frames','color',h.xlabelcolor)
        end
        
%         for ii=1:size(x,1)
%                 D(ii,:)=sqrt( (x(ii)-x(:)).^2 +(y(ii)-y(:)).^2 +(z(ii)-z(:)).^2 );
%         end
%         
%         grain=1;
%         xmin=min(x);ymin=min(y);zmin=min(z);
%         I3=zeros(grain*ceil(max(x)-xmin+10),grain*ceil(max(y)-ymin+10),ceil(max(z)-min(z))+10);
% %         D2=pdist2(x,y);
%         for ii=1:size(x,1)
%             I3( 5+round((x(ii)-xmin)), 5+round((y(ii)-ymin)),5+round((z(ii)-zmin)))=c(ii);
%         end
        if get(h.chb_voxels,'Value')==1
            makeVol_vPALM(x,y,z,str2num(get(h.Nvoxels,'String')),0,0,h);
        else
            scatter3(x,y,z,4,c)
            xlabel('X,px','color',h.xlabelcolor)
            ylabel('Y,px','color',h.xlabelcolor)
            zlabel('Intensity, a.u.','color',h.xlabelcolor)
            colorbar
            set(h.axes2,'Xcolor',[0.5 0.5 0.5]);
            set(h.axes2,'Ycolor',[0.5 0.5 0.5]);
            view(str2num(get(h.AZangle,'String')),str2num(get(h.ELangle,'String')) );
        end
    end
    
    status(h,'g',strcat('2d data plotted in :',num2str(toc,2),' s'));

%     colormap parameters.colormap;

%% density2d
function density=density2d(h,x,y)
status(h,'r','working...')
% get(h.h_sliders.Xmax_slider,'Max')
% det_thresh=str2num(get(h.detect_thresh,'String'))./str2num(get(h.pixelsize,'String'));
det_thresh=str2num(get(h.h_sliders.detect_thresh,'String'))./str2num(get(h.pixelsize,'String'));

parfor iparticle=1:size(x,1)
    n= find( sqrt( (x(iparticle)-x).^2 + (y(iparticle)-y).^2 ) < det_thresh);
    density(iparticle)=size(n,1 );
    
end

disp('density calculated');
    status(h,'g','ok, i am done')

%% density3d
function density=density3d(h,x,y,z)
status(h,'r','working...');
tic
% det_thresh=str2num(get(h.detect_thresh,'String'))./str2num(get(h.pixelsize,'String'));
det_thresh=str2num(get(h.h_sliders.detect_thresh,'String'))./str2num(get(h.pixelsize,'String'));

parfor iparticle=1:size(x,1)
    n=find( sqrt( (x(iparticle)-x).^2 + (y(iparticle)-y).^2 + (z(iparticle)-z).^2) < det_thresh) ;
    density(iparticle)=size(n,1);
    
end

disp('density calculated');
status(h,'g',strcat('ok, i am done in ',num2str(toc),' s'));



%% reconstruction
function plot_reconstruction(hObject, eventdata, h)
status(h,'g','Calculating reconstruction...');
m=h.m;
pointing_precision_px=str2num(get(h.h_sliders.detect_thresh,'String'))./str2num(get(h.pixelsize,'String'));

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

    end
end

if get(h.chb_DIPimage,'Value')==0
    px = str2num(get(h.pixelsize,'String'));
    x = x*px;
    y = y*px;
    uniform_peaks=get(h.chb_uniform_peaks,'Value');
    border=10;
    Pixel_to_photon=1;
    max_intensity=1;
    hot=1;
    Contrast_tool=1;
    ps=str2num(get(h.stepsize_reconstruction,'String'));
    % [I,xcoor,ycoor,Imax]=vPALM_reconstruction_v5(x,y,s,A, ps,10,1,1,1,1,logyes);
    [PALM_image]=vPALM2_analysis_image_reconstruction_v4b(x, y, s, A, ps, border, Pixel_to_photon,uniform_peaks, max_intensity, hot, Contrast_tool);
else
    % DIPimage
    ps=str2num(get(h.stepsize_reconstruction,'String'));
    px = str2num(get(h.pixelsize,'String'));
    
%     coords(:,1) = x-min(x);
%     coords(:,2) = y-min(y);
%     coords(:,3) = t;
%     uniform_peaks=get(h.chb_uniform_peaks,'Value');
    s=str2num(get(h.resolution_reconstruction,'String'))./str2num(get(h.pixelsize,'String'));

    vPALM2_reconstruction_time(x,y,t,s,ps,px)
%     superzoom = ceil(px/ps);
%     szx = superzoom * ceil(max(coords(:,1)));
%     szy = superzoom * ceil(max(coords(:,2)));
% 
%     im = binlocalizations(coords, szx, szy, superzoom);
%     h1=dipshow(im);
%     dipmapping(h1,[0 5]);%,'colormap',hot);


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


%% load_MTT_ch1
% function load_MTT_ch1(hObject, eventdata, h)
% 
% axes(h.axes1);
% cla;
% 
% status(h,'c','Give the filename of the image file');
% [h.fullFileName,path1] = uigetfile('*.mat','Open MTT File');
% cd(path1);
% load(h.fullFileName); % reads calfit and phase_fromto
% 
% set(h.filename,'String',h.fullFileName);
% 
% h.m=matrice_results;
% % Update handles structure
% setcallbacks(h)
% 
% %% load_MTT_ch2
% function load_MTT_ch2(hObject, eventdata, h)
% 
% status(h,'c','Give the filename of the image file');
% [h.fullFileName,path1] = uigetfile('*.mat','Open MTT File');
% cd(path1);
% load(h.fullFileName); % reads calfit and phase_fromto
% 
% set(h.filename,'String',h.fullFileName);
% 
% h.m2=matrice_results;
% 
% % Update handles structure
% setcallbacks(h)
% 
% %% plot_localizations_ch2
% function plot_localizations_ch1(hObject, eventdata, h)
% 
% m=h.m;
% 
% % axes(handles.axes1);
% % cla;
% 
% figure(1)
% hold on
% 
% plot(m(2,:),m(3,:),'or','MarkerSize',12)
% grid
% 
% 
% %% load_MTT_ch2
% function plot_localizations_ch2(hObject, eventdata, h)
% 
% m2=h.m2;
% 
% figure(1)
% hold on
% plot(m2(2,:),m2(3,:),'+g')
% grid


%% density image
function plot_reconstruction(hObject, eventdata, h)
status(h,'g','Calculating density...');
m=h.m;

pointing_precision_px=str2num(get(h.h_sliders.detect_thresh,'String'))./str2num(get(h.pixelsize,'String'));




%% 
function m2=filter_sliders(h,m)
    p=read_parameters(h);
    
    if h.parameters.fit_3dxcorr==0 & h.parameters.fit_3dgaussian==0
        
        if get(h.chb_drift_correction,'Value')==0

            r=find(m(:,1)>p.Frmin & m(:,1)<p.Frmax & m(:,2)>p.Xmin & m(:,2)<p.Xmax ...
                    & m(:,3)>p.Ymin & m(:,3)<p.Ymax & m(:,4)>p.Imin & m(:,4)<p.Imax);
            m2=m(r,:);   
        else %uses drift correction
            m=h.m_drift;
            r=find(m(:,1)>p.Frmin & m(:,1)<p.Frmax & m(:,2)>p.Xmin & m(:,2)<p.Xmax ...
                    & m(:,3)>p.Ymin & m(:,3)<p.Ymax & m(:,4)>p.Imin & m(:,4)<p.Imax);
            m2=m(r,:);   

        end
%         r1=find(m(:,1)>p.Frmin & m(:,1)<p.Frmax);
%         r2=find(m(r1,4)>p.Imin & m(r1,4)<p.Imax);
%         r3=find(m(r2,2)>p.Xmin & m(r2,2)<p.Xmax);
%         r4=find(m(r3,3)>p.Ymin & m(r3,3)<p.Ymax);
%         m2=m(r4,:);


    elseif h.parameters.fit_3dxcorr==1 %      plots 3d xcorr localizations

        r=find(m(:,1)>p.Frmin & m(:,1)<p.Frmax & m(:,2)>p.Xmin & m(:,2)<p.Xmax ...
            & m(:,3)>p.Ymin & m(:,3)<p.Ymax & m(:,4)>p.Imin & m(:,4)<p.Imax ...
            & m(:,13)>p.Zmin & m(:,13)<p.Zmax);
        m2=m(r,:);

    elseif h.parameters.fit_3dgaussian==1 %      plots 3d gaussian localizations    
   
        if h.parameters.v_gauss_cal==0 % no calibation file, wil show ratio wy/wx;
            z=m(:,13)./m(:,14);
            m(:,15)=z;
            disp('plotting wx/wy as no cal is available')
        else % uses cal
%             z=calc_z_from_cal(h,m,h.parameters);
%             m(:,15)=z;
            if get(h.chb_3dcalapplied,'Value')==1
                m=h.m_3d;
                disp('using calibration')
            else
                m=h.m;
                z=m(:,13)./m(:,14);
                m(:,15)=z;
                disp('no calibration') 
            end
        end

        r=find(m(:,1)>p.Frmin & m(:,1)<p.Frmax & m(:,2)>p.Xmin & m(:,2)<p.Xmax ...
            & m(:,3)>p.Ymin & m(:,3)<p.Ymax & m(:,4)>p.Imin & m(:,4)<p.Imax ...
            & m(:,15)>p.Zmin & m(:,15)<p.Zmax);
        m2=m(r,:);

end



function offsetx_slider(hObject, eventdata, h)

%     set(h.offsetx_slider,'String',num2str(get(h.offsetx_slider,'Value')));
    plot_localizations_2d(hObject,eventdata,h)

function offsety_slider(hObject, eventdata, h)

%     set(h.detect_thresh,'String',num2str(get(h.detect_thresh_slider,'Value')));
    plot_localizations_2d(hObject,eventdata,h)

    
    
function detect_thresh(hObject, eventdata, h)
    set(h.detect_thresh_slider,'Value',str2num(get(h.detect_thresh,'String')));
    plot_localizations_2d(hObject,eventdata,h)


%%
function Zcolormap(hObject, eventdata, h)

v=get(h.Zcolormap,'Value');

switch v
    case 1
        colormap jet
    case 2
        colormap hot
    case 3
        colormap HSV
    case 4
        colormap gray
end

% plot_localizations(hObject,eventdata,h);

%%
function zplotoption(hObject, eventdata, h)

v=get(h.zplotoption,'Value');

plot_localizations_3d(hObject,eventdata,h)


%% make_roi
function make_roi(hObject, eventdata, h)


axes(h.axes1)
h1 = imrect(h.axes1);
rectzoom = wait(h1); % wait until you double-click on the image and return the position of the rectangle

h.rectzoom=rectzoom;

set(h.h_sliders.Xmax,'String',num2str(h.rectzoom(1)+rectzoom(3)));
set(h.h_sliders.Xmin,'String',num2str(h.rectzoom(1)));
set(h.h_sliders.Ymax,'String',num2str(h.rectzoom(2)+rectzoom(4)));
set(h.h_sliders.Ymin,'String',num2str(h.rectzoom(2)));

if h.rectzoom(1)+rectzoom(3)>get(h.h_sliders.Xmax_slider,'Max')
    set(h.h_sliders.Xmax_slider,'Max',h.rectzoom(1)+rectzoom(3)+1)
end
set(h.h_sliders.Xmax_slider,'Value',h.rectzoom(1)+rectzoom(3));

if h.rectzoom(1)<get(h.h_sliders.Xmin_slider,'Min')
    set(h.h_sliders.Xmin_slider,'Min',h.rectzoom(1)-1)
end
set(h.h_sliders.Xmin_slider,'Value',h.rectzoom(1));

if h.rectzoom(2)+rectzoom(4)>get(h.h_sliders.Ymax_slider,'Max')
    set(h.h_sliders.Ymax_slider,'Max',h.rectzoom(2)+rectzoom(4)+1)
end
set(h.h_sliders.Ymax_slider,'Value',h.rectzoom(2)+rectzoom(4));

if h.rectzoom(2)<get(h.h_sliders.Ymin_slider,'Min')
    set(h.h_sliders.Ymin_slider,'Min',h.rectzoom(2)-1)
end
set(h.h_sliders.Ymin_slider,'Value',h.rectzoom(2));

setcallbacks(h);


plot_localizations_2d(hObject,eventdata,h)


%% reset_roi
function reset_roi(hObject, eventdata, h)

slider_setup(h,h.m,h.parameters)
plot_localizations_2d(hObject,eventdata,h)

%% save_axes1
function save_axes1(hObject, eventdata, h)

axes(h.axes1)
% export_fig test3.png
F=getframe(h.axes1);               %select axes in GUI
figure();                                          %new figure
image(F.cdata);  %show selected axes in new figure
saveas(gcf, strcat(h.fullFileName,'_axes1.fig'), 'fig');                    %save figure
close(gcf); %and close it

% 
% axes(h.axes1)
% h1=gcf;`
% h2=figure;
% objects=allchild(h1);
% copyobj(get(h1,'children'),h2);

%% save_axes2
% function save_axes2(hObject, eventdata, h)
% 
% axes(h.axes2)
% % export_fig test3.png
% F=getframe(h.axes2);               %select axes in GUI
% figure();                                          %new figure
% image(F.cdata);  %show selected axes in new figure
% saveas(gcf, strcat(h.fullFileName,'_axes2.png'), 'png');                    %save figure
% close(gcf);                                       %and close it

%% selectbeads
function selectbeads(hObject, eventdata, h)

h=vPALM2_Beads_detection(h)
m2=h.m;
m2(:,2)=m2(:,2)-h.Ref_position(m2(:,1),1);
m2(:,3)=m2(:,3)-h.Ref_position(m2(:,1),2);
h.m_drift=m2;

set(h.chb_drift_correction,'Value',1);
setcallbacks(h);
disp('done');

% make_roi(hObject, eventdata, h);
% m=filter_sliders(h,h.m);
% 
% axes(h.axes2);
% 
% I=m(:,4);
% maxI=max(I);
% 
% select=find(I>0);%get(h.showbeads_slider,'Value')*maxI);
% t=m(select,1);
% x=m(select,2);
% y=m(select,3);
% meanx=mean(x);
% meany=mean(y);
% 
% xm=meanfilter(100,x)'+meanx;
% ym=meanfilter(100,y)'+meany;
% 
% 
% 
% %makes interpolation
% frmin=min(h.m(:,1));
% frmax=max(h.m(:,1));
% frrange=[frmin:1:frmax+1];
% h.tm=frrange;
% xm=interp1(t,xm,frrange);
% ym=interp1(t,ym,frrange);
% h.xm(frrange)=xm';
% h.ym(frrange)=ym';
% 
% x2=x-h.xm(t)';
% y2=y-h.ym(t)';
% 
% % for i=1:size(x,1)
% %     x2(:)=x(:)-h.xm(find(h.tm==t(:)));
% %     y2(:)=y(:)-h.ym(find(h.tm==t(:)));
% % end
% 
% hold off,
% plot(x2,y2,'.b')
% 
% x2=x-h.xm(t)';
% y2=y-h.ym(t)';
function apply_3dcal(hObject, eventdata, h)

if get(h.chb_cal_3d,'Value')==1
    h.m_3d=PALM2_calc_3d(h);

    set(h.chb_3dcalapplied,'Value',1);
else
    disp('No calibration was done!');
end

% Update handles structure
setcallbacks(h)

% %% showbeads
% function showbeads(hObject, eventdata, h)
% 
% m=filter_sliders(h,h.m);
% plot_localizations_2d(hObject,eventdata,h)
% 
% axes(h.axes1);
% hold on,
% I=m(:,4);
% 
% maxI=max(I)
% select=find(I>get(h.showbeads_slider,'Value')*maxI);
% 
% x=m(select,2);
% y=m(select,3);
% plot(x,y,'*b','MarkerSize',8)

% still to program...


%% saves data to file
function save_data(hObject, eventdata, h)

disp_nice('Give the filename of the mat file')
% write_parameters(hObject,h,parameters);

filename=strcat(h.fullFileName,'_vPALM.mat');

switch get(h.load_options,'Value')
    
    case 2 %qPALM
%         [fullFileName,path1] = uigetfile('*_analysed.mat','Open MAT File');
%         cd(path1);
%         load(fullFileName);
%         h.m=matrice_results';
%         h.fullFileName=parameters.fullFileName;
%         h.parameters=parameters;
%         
%         if parameters.v_gauss_cal==1
%           set(h.chb_cal_3d,'Value',1);
%         else
%           set(h.chb_cal_3d,'Value',0);
%         end
    case 3 % PALMcbs
%         [fullFileName,path1] = uigetfile('*_analysed.mat','Open MAT File');
%         cd(path1);
%         load(fullFileName);
%         h.m=AllData.m1';
%         h.m(:,2)=AllData.m1(3,:);
%         h.m(:,3)=AllData.m1(2,:);
% 
%         h.fullFileName=fullFileName;
%         parameters.fit_2dgaussian=1;
%         parameters.fit_3dgaussian=0;
%         parameters.fit_3dxcorr=0;
%         h.parameters=parameters;
%         if get(h.chb_image_load,'Value')==1
%             disp_nice('Give the filename of the image file')
%             [TIFfilename,path1] = uigetfile('*.tif','Open TIF File');
%             h.fullFileName=TIFfilename;
%             parameters.fullFileName=TIFfilename;
%         else
%             h.fullFileName=fullFileName; 
%             parameters.fullFileName=fullFileName;
%         end

    
    case 4 %MTT
%         [fullFileName,path1] = uigetfile('*_analysed.mat','Open MAT File');
%         cd(path1);
%         load(fullFileName);
%         h.m=matrice_results';    
%         parameters.fit_2dgaussian=1;
%         parameters.fit_3dgaussian=0;
%         parameters.fit_3dxcorr=0;
%         
%         h.parameters=parameters;
%         if get(h.chb_image_load,'Value')==1
%             disp_nice('Give the filename of the image file')
%             [TIFfilename,path1] = uigetfile('*.tif','Open TIF File');
%             h.fullFileName=TIFfilename;
%             parameters.fullFileName=TIFfilename;
%         else
%             h.fullFileName=fullFileName; 
%             parameters.fullFileName=fullFileName;
%         end

    case 5 %ENS
%         [fullFileName,path1] = uigetfile('*_analysed.mat','Open MAT File');
%         cd(path1);
%         load(fullFileName);
%         h.m=matrice_results';
%         h.fullFileName=fullFileName;
%         parameters.fit_2dgaussian=0;
%         parameters.fit_3dgaussian=0;
%         parameters.fit_3dxcorr=1;
%         h.m(:,7:14)=[zeros(size(h.m,1),6),h.m(:,4),zeros(size(h.m,1),1)];
%         h.m(:,4)=h.m(:,5);
%         parameters.fullFileName=fullFileName;
%         h.parameters=parameters;
    case 6 %micromanager
        
        if get(h.chb_drift_correction,'Value')==0
            m0=h.m;
        else
            m0=h.m_drift;
        end
            
        
        m=filter_sliders(h,m0);

        if get(h.chb_saveappend,'Value')==1
            clear data;
            load(filename);
            data.Nframe=[data.Nframe; m(:,1)]; % frame number
            data.x=[data.x;(m(:,2)-1).*str2num(get(h.pixelsize,'String'))]; % pos x
            data.y=[data.y;(m(:,3)-1).*str2num(get(h.pixelsize,'String'))]; % pos y
            data.Intensity=[data.Intensity;m(:,4)]; % intensity
            data.Resolution=[data.Resolution;m(:,5)]; % resolution, nm
            data.Backround=[data.Backround;m(:,6)]; %background
            data.Molecule_id=[data.Molecule_id;m(:,7)]; %molecule
            data.Molecule2=[data.Molecule2;m(:,8)]; %molecule
            data.Channel=[data.Channel;m(:,9)]; %channel
            data.Slice=[data.Slice;m(:,10)]; %slice
            data.Width=[data.Width;m(:,11)]; %width
            data.A= [data.A;m(:,12)]; %a ???
            data.Theta= [data.Theta;m(:,13)]; %theta
            data.Xpos= [data.Xpos;m(:,14)]; %xpos
            data.Ypos=[data.Ypos;m(:,15)]; %ypos
        else
            clear data;
            data.Nframe=m(:,1); % frame number
            data.x=(m(:,2)-1).*str2num(get(h.pixelsize,'String')); % pos x
            data.y=(m(:,3)-1).*str2num(get(h.pixelsize,'String')); % pos y
            data.Intensity=m(:,4); % intensity
            data.Resolution=m(:,5); % resolution, nm
            data.Backround=m(:,6); %background
            data.Molecule_id=m(:,7); %molecule
            data.Molecule2=m(:,8); %molecule
            data.Channel=m(:,9); %channel
            data.Slice=m(:,10); %slice
            data.Width=m(:,11); %width
            data.A= m(:,12); %a ???
            data.Theta= m(:,13); %theta
            data.Xpos= m(:,14); %xpos
            data.Ypos=m(:,15); %ypos
        end
      
        clear parameters;
        parameters.fullFileName=h.fullFileName;
        parameters.fit_2dgaussian=h.parameters.fit_2dgaussian;
        parameters.fit_3dgaussian=h.parameters.fit_3dgaussian;
        parameters.fit_3dxcorr=h.parameters.fit_3dxcorr;
        parameters.fullFileName=h.parameters.fullFileName;
        parameters.pixelsize=str2num(get(h.pixelsize,'String'));
        parameters.chb_drift_correction=get(h.chb_drift_correction,'Value');
        parameters.chb_image_load=get(h.chb_image_load,'Value');
        parameters.chb_cal_3d=get(h.chb_cal_3d,'Value');
        parameters.chb_3dcalapplied=get(h.chb_3dcalapplied,'Value');
        parameters.v_gauss_cal=h.parameters.v_gauss_cal;

   
        save(filename,'data','parameters');

        
 case 7 %micromanager

 if get(h.chb_drift_correction,'Value')==0
            m0=h.m;
        else
            m0=h.m_drift;
        end
            
        if get(h.chb_cal_3d,'Value')==0
            m=filter_sliders(h,m0);
        else
            m=filter_sliders(h,h.m_3d);       
        end
        
        if get(h.chb_saveappend,'Value')==1
            clear data;
            load(filename);
            data.Nframe=[data.Nframe; m(:,1)]; % frame number
            data.x=[data.x;(m(:,2)-1).*str2num(get(h.pixelsize,'String'))]; % pos x
            data.y=[data.y;(m(:,3)-1).*str2num(get(h.pixelsize,'String'))]; % pos y
            data.Intensity=[data.Intensity;m(:,4)]; % intensity
            data.Resolution=[data.Resolution;m(:,5)]; % resolution, nm
            data.Backround=[data.Backround;m(:,6)]; %background
            data.Molecule_id=[data.Molecule_id;m(:,7)]; %molecule
            data.Molecule2=[data.Molecule2;m(:,8)]; %molecule
            data.Channel=[data.Channel;m(:,9)]; %channel
            data.Slice=[data.Slice;m(:,10)]; %slice
            data.Width=[data.Width;m(:,11)]; %width
            data.A= [data.A;m(:,12)]; %a ???
            data.Theta= [data.Theta;m(:,13)]; %theta
            data.Xpos= [data.Xpos;m(:,14)]; %xpos
            data.Ypos=[data.Ypos;m(:,15)]; %ypos
            data.Wx=[data.Wx;m(:,13)]; %ypos
            data.Wy=[data.Wy;m(:,14)]; %ypos
            data.Asymmetry=[data.Asymmetry;m(:,11)]; %ypos
            data.z=[data.z;m(:,15)]; %ypos   

        else
            clear data;
            data.Nframe=m(:,1); % frame number
            data.x=(m(:,2)-1).*str2num(get(h.pixelsize,'String')); % pos x
            data.y=(m(:,3)-1).*str2num(get(h.pixelsize,'String')); % pos y
            data.Intensity=m(:,4); % intensity
            data.Resolution=m(:,5); % resolution, nm
            data.Backround=m(:,6); %background
            data.Molecule_id=m(:,7); %molecule
            data.Molecule2=m(:,8); %molecule
            data.Channel=m(:,9); %channel
            data.Slice=m(:,10); %slice
            data.Width=m(:,11); %width
            data.A= m(:,12); %a ???
            data.Theta= m(:,13); %theta
            data.Xpos= m(:,14); %xpos
            data.Ypos=m(:,15); %ypos
            data.Wx=m(:,13); %ypos
            data.Wy=m(:,14); %ypos
            data.Asymmetry=m(:,11); %ypos
            data.z=m(:,15); %ypos            
        end
      
        clear parameters;
        parameters.fullFileName=h.fullFileName;
        parameters.fit_2dgaussian=h.parameters.fit_2dgaussian;
        parameters.fit_3dgaussian=h.parameters.fit_3dgaussian;
        parameters.fit_3dxcorr=h.parameters.fit_3dxcorr;
        parameters.fullFileName=h.parameters.fullFileName;
        parameters.pixelsize=str2num(get(h.pixelsize,'String'));
        parameters.chb_drift_correction=get(h.chb_drift_correction,'Value');
        parameters.chb_image_load=get(h.chb_image_load,'Value');
        parameters.chb_cal_3d=get(h.chb_cal_3d,'Value');
        parameters.chb_3dcalapplied=get(h.chb_3dcalapplied,'Value');
        parameters.v_gauss_cal=h.parameters.v_gauss_cal;
        
        if get(h.chb_cal_3d,'Value')==1
            parameters.wxPoly = h.wxPoly ;
            parameters.wyPoly = h.wyPoly ;
            parameters.ccPolyWtoZ = h.ccPolyWtoZ ;
        end
        
        save(filename,'data','parameters');

end

        
status(h,'g','everything went apparently fine')
disp(strcat('Data saved to: ',filename))

function calc_res(hObject, eventdata, h)

status(h,'g','Calculating image resolution')
% disp(strcat('Data saved to: ',filename))

h.m2=filter_sliders(h,h.m);


ps=str2num(get(h.stepsize_reconstruction,'String'));
px = str2num(get(h.pixelsize,'String'));

x=h.m2(:,2);
y=h.m2(:,3);
t=h.m2(:,1);

coords(:,1) = x-min(x);
coords(:,2) = y-min(y);
coords(:,3) = t;

superzoom = ceil(px/ps);
szx = superzoom * ceil(max(coords(:,1)));
szy = superzoom * ceil(max(coords(:,2)));

im = binlocalizations(coords, szx, szy, superzoom);
hh=dipshow(im);
dipmapping(hh,[0 5],'colormap',hot)

fprintf('\n -- computing FIRE --\n')
[fire_value, ~, fireH, fireL] = postoresolution(coords, szx, superzoom); % in super-resolution pixels
fprintf('FIRE value %2.1f +- %2.2f [px]\n', fire_value, (fireL-fireH)/2);
fprintf('FIRE value %2.1f +- %2.2f [nm]\n', fire_value*px/superzoom, (fireL-fireH)/2*px/superzoom);

status(h,'g',strcat('FIRE value :',num2str(fire_value*px/superzoom),' +- ',num2str((fireL-fireH)/2*px/superzoom))) 


function save_PALMvis(hObject, eventdata, h)

vPALM2_save_to_PALMvis(h);




function stats(hObject, eventdata, h)

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
        status(h,'g','calculating N local/px statistics');
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
        res=m(:,4);
        hist((res))
        disp('Intensity');
        xlabel('Intensity','color',h.xlabelcolor)
        ylabel('PDF (a.u.)','color',h.xlabelcolor)
        grid
        set(gca,'Xcolor',[0.5 0.5 0.5]);
        set(gca,'Ycolor',[0.5 0.5 0.5]);
        status(h,'g',strcat('mean res',num2str(mean(res)),' nm'));        
               
        
end



%%  Utility functions for MYGUI          

%% read_parameters
function p=read_parameters(h)
    p.Zmax=str2num(get(h.h_sliders.Zmax,'String'));
    p.Zmin=str2num(get(h.h_sliders.Zmin,'String'));
    p.Xmax=str2num(get(h.h_sliders.Xmax,'String'));
    p.Xmin=str2num(get(h.h_sliders.Xmin,'String'));
    p.Ymax=str2num(get(h.h_sliders.Ymax,'String'));
    p.Ymin=str2num(get(h.h_sliders.Ymin,'String'));
    p.Frmax=str2num(get(h.h_sliders.Frmax,'String'));
    p.Frmin=str2num(get(h.h_sliders.Frmin,'String'));
    p.Imax=str2num(get(h.h_sliders.Imax,'String'));
    p.Imin=str2num(get(h.h_sliders.Imin,'String'));
    
%     set(h.h_sliders.detect_thresh,'String','-10');



%% set callbacks    
function setcallbacks(h) 

set(h.push_load_file,'callback',{@load_file,h});
set(h.plot_localizations_3d,'callback',{@plot_localizations_3d,h});
set(h.plot_localizations_2d,'callback',{@plot_localizations_2d,h});

set(h.plot_reconstruction,'callback',{@plot_reconstruction,h});
set(h.plot_density,'callback',{@plot_density,h});


set(h.cal_3d,'callback',{@cal_3d,h});

set(h.calc_res,'callback',{@calc_res,h});
set(h.Zcolormap,'callback',{@Zcolormap,h});
set(h.zplotoption,'callback',{@zplotoption,h});
set(h.stats,'callback',{@stats,h});

set(h.save_data,'callback',{@save_data,h});


set(h.apply_3dcal,'callback',{@apply_3dcal,h});
set(h.save_PALMvis,'callback',{@save_PALMvis,h});

set(h.make_roi,'callback',{@make_roi,h});
set(h.reset_roi,'callback',{@reset_roi,h});
set(h.save_axes1,'callback',{@save_axes1,h});
set(h.offsetx_slider,'callback',{@offsetx_slider,h});
set(h.offsety_slider,'callback',{@offsety_slider,h});
set(h.selectbeads,'callback',{@selectbeads,h});
% set(h.showbeads,'callback',{@showbeads,h});



%%

function m=filter_3dgauss(m0,parameters,h)

if get(h.load_options,'Value')==6

    ellip=m0(:,11);
    
    mean_w=m0(:,13)/h.pixelsize;

else
    
    ellip=m0(:,13)./m0(:,14);
    
    mean_w=sqrt(m0(:,13).*m0(:,14));
    
end

    keep_values=find(ellip>parameters.min_ellip & ellip<parameters.max_ellip & mean_w>parameters.min_w & mean_w< parameters.max_w);

    m=m0(keep_values,:);
  
function status(h,color,text)
disp(text);
set(h.status,'String',text,'ForegroundColor',color);


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
