%% load_file
function h=vPALM_load_file(hObject,eventdata,h)

disp('Give the filename of the mat file')
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
            disp('Give the filename of the image file')
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
            disp('Give the filename of the image file')
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
        A=importdata(fullFileName);
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
        
    case 'ThunderSTORM' %ThunderSTORM
        [fullFileName,path1] = uigetfile('*.*','Open ThunderSTORM File');
        cd(path1);
        A=importdata(fullFileName);
        
        clear_m_matrices;
        %1 if
        %2 frame number
        %3,4 x and y in nm
        %5 sigma [nm]
        %6 intensity [photon]
        %7 offset [photon]
        %8 bkgstd [photon]
        %9 chi2
        %10 uncertainty [nm]
        if size(A.data,2)==10 % 2D data
            m(:,1)=A.data(:,2); % frame number
            m(:,2)=A.data(:,3)./str2num(get(h.pixelsize,'String')); % pos x in px
            m(:,3)=A.data(:,4)./str2num(get(h.pixelsize,'String')); %pos y in px
            m(:,4)=A.data(:,6); % intensity
            m(:,5)=A.data(:,10); % resolution, nm
            m(:,6)=A.data(:,8); % background
            m(:,7)=A.data(:,9); % fit residue, chi square
            m(:,8)=A.data(:,5); %sigma1 in nm
            
            m(:,13)=zeros(size(A.data(:,1))); %Wx %%added tue1405
            m(:,14)=zeros(size(A.data(:,1))); %Wy %%added tue1405
            m(:,10)=zeros(size(A.data(:,1))); %angle
            m(:,11)=zeros(size(A.data(:,1))); %asymmetry
            m(:,12)=zeros(size(A.data(:,1))); %channel
            m(:,9)=zeros(size(A.data(:,1))); %pos

        elseif  size(A.data,2)==12 % 2D data
            m(:,1)=A.data(:,2); % frame number
            m(:,2)=A.data(:,3)./str2num(get(h.pixelsize,'String')); % pos x in px
            m(:,3)=A.data(:,4)./str2num(get(h.pixelsize,'String')); %pos y in px
            m(:,4)=A.data(:,8); % intensity
            m(:,5)=A.data(:,12); % resolution, nm
            m(:,6)=A.data(:,10); %background
            m(:,7)=A.data(:,11); % fit residue, chi square
            m(:,8)=A.data(:,6); %sigma1 in nm
            m(:,9)=A.data(:,7); %sigma2 in nm
            m(:,10)=zeros(size(A.data(:,1))); %angle
            m(:,11)=zeros(size(A.data(:,1))); %asymmetry
            m(:,12)=zeros(size(A.data(:,1))); %channel
            m(:,13)=zeros(size(A.data(:,1))); %Wx %%added tue1405
            m(:,14)=zeros(size(A.data(:,1))); %Wy %%added tue1405
            m(:,15)=A.data(:,5); %z in nm from calibration in ThunderSTORM
        end
          
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
        
        if size(A.data,2)==10 % 2D data
            h.fullFileName=fullFileName;
            h.parameters.fit_2dgaussian=1;
            h.parameters.fit_3dgaussian=0;
            h.parameters.fit_3dxcorr=0;
        elseif size(A.data,2)==12 % 3D data        
            h.fullFileName=fullFileName;
            h.parameters.fit_2dgaussian=0;
            h.parameters.fit_3dgaussian=1;
            h.parameters.fit_3dxcorr=0;
            h.parameters.v_gauss_cal=1;
            h.parameters.min_ellip=0;
            h.parameters.max_ellip=MAX_ELLIP;
            h.parameters.min_w=0 ;
            h.parameters.max_w=MAX_r0*str2num(get(h.pixelsize,'String'));
            set(h.chb_cal_3d,'Value',1);
            set(h.chb_3dcalapplied,'Value',1);
            h.m_3d=m2;

        end

        
        h.parameters.fullFileName=fullFileName;
%         h.parameters=parameters;
        set(h.chb_image_load,'Value',0);
end

        
status(h,'g','everything went apparently fine')

set(h.filename,'String',fullFileName,'ForegroundColor',[0 0 1]);
set(h.chb_drift_correction,'Value',0);

set(h.Nparticles_txt,'String',num2str(size(h.m,1)));

set(h.chb_data_loaded,'Value',1);

slider_setup(h,h.m, h.parameters);

       
% Update handles structure
setcallbacks(h)
