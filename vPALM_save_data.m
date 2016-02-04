%% saves data to file
function vPALM_save_data(hObject, eventdata, h)

disp('Give the filename of the mat file')
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
%         
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
            data.Backround=m(:,5); %background
            data.Molecule_id=m(:,5); %molecule
            data.Molecule2=m(:,5); %molecule
            data.Channel=m(:,5); %channel
            data.Slice=m(:,5); %slice
            data.Width=m(:,5); %width
            data.A= m(:,5); %a ???
            data.Theta= m(:,5); %theta
            data.Xpos= m(:,5); %xpos
            data.Ypos=m(:,5); %ypos
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
%         parameters.v_gauss_cal=h.parameters.v_gauss_cal;

   
        save(filename,'data','parameters');

        

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
    case 6 %micromanager 2d
        
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

        
 case 7 %micromanager 3d

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
