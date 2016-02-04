%% 2d-3d commands for figure3D


function h_vis = vPALM2_2d3d_visualization(h,varargin)

disp('ass')

h = get(h.MainFig, 'UserData');  
h_vis.h=h;

% h=varargin{1}; % handles h
% h_vis.axes2=h.axes2;
scrsz = get(0,'ScreenSize');
% FIGH=scrsz(4)/1;
% FIGW=scrsz(3)*.8;
FIGH=250;
FIGW=300; %1092

h_vis.axes2=h.axes2; %figure();

% axis2_size=[1 1 250 250];
% 
% h_vis.axes2= axes('Parent',h_vis.Figure3D,'Position',axis2_size);
% 
MainWindowPos=get(h.MainFig,'Position')

h_vis.fig = figure(...
              'Name','2d/3d vis',...
              'Units','pixels',...
              'MenuBar','none',...
              'Toolbar','none',...
              'Color','w',...              
              'Position',[MainWindowPos(3)+5 MainWindowPos(4)-FIGH+100 FIGW FIGH],...
              'Visible','on',...
              'CloseRequestFcn',@newclosefunction);

get(h_vis.fig,'Position')

%%  Construct the components
h_vis=vPALM2_2d3d_visualization_gui(h_vis,FIGH,FIGW);

vPALM2_2d3d_visualization_callbacks(h_vis);

end


%% callbacks
function vPALM2_2d3d_visualization_callbacks(h_vis)

set(h_vis.xyplane,'callback',{@xyplane,h_vis});
set(h_vis.yzplane,'callback',{@yzplane,h_vis});
set(h_vis.xzplane,'callback',{@xzplane,h_vis});
set(h_vis.xyzplane,'callback',{@xyzplane,h_vis});
set(h_vis.rotAZ,'callback',{@rotAZ,h_vis});
set(h_vis.rotEL,'callback',{@rotEL,h_vis});
set(h_vis.rotAZ_neg,'callback',{@rotAZ_neg,h_vis});
set(h_vis.rotEL_neg,'callback',{@rotEL_neg,h_vis});


end

%% functions
function newclosefunction(hObject, eventdata)
disp('Remember to close this GUI using the main vPALM window!');
end

function xyplane(hObject, eventdata, h_vis)
axes(h_vis.axes2);
view(2);
end

function yzplane(hObject, eventdata, h_vis)
axes(h_vis.axes2);
view([1 0 0]);
end

function xzplane(hObject, eventdata, h_vis)
axes(h_vis.axes2);
view([0 1 0]);
end

function xyzplane(hObject, eventdata, h_vis)
axes(h_vis.axes2);
view(3);
end


function rotAZ(hObject, eventdata, h_vis)

h = get(h_vis.h.MainFig, 'UserData');  
set(h.AZangle,'String',num2str( str2num(get(h.AZangle,'String')) +10 ));
axes(h_vis.axes2);
view(str2num(get(h.AZangle,'String')),str2num(get(h.ELangle,'String')) );
end

function rotEL(hObject, eventdata, h_vis)
h = get(h_vis.h.MainFig, 'UserData');  

set(h.ELangle,'String',num2str( str2num(get(h.ELangle,'String')) +10 ));
axes(h_vis.axes2);
view(str2num(get(h.AZangle,'String')),str2num(get(h.ELangle,'String')) );
end

function rotAZ_neg(hObject, eventdata, h_vis)
h = get(h_vis.h.MainFig, 'UserData');  

set(h.AZangle,'String',num2str( str2num(get(h.AZangle,'String')) -10 ));
axes(h_vis.axes2);
view(str2num(get(h.AZangle,'String')),str2num(get(h.ELangle,'String')) );
end

function rotEL_neg(hObject, eventdata, h_vis)
h = get(h_vis.h.MainFig, 'UserData');  

set(h.ELangle,'String',num2str( str2num(get(h.ELangle,'String')) -10 ));
axes(h_vis.axes2);
view(str2num(get(h.AZangle,'String')),str2num(get(h.ELangle,'String')) );
end

