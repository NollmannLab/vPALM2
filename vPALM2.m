% vPALM GUI
function varargout = vPALM2(varargin)

%% M Nollmann, Sep 2015.
% Centre de Biochimie Structurale
% CNRS, 29 rue de Navacelles, Montpellier, France

% vPALM Brief description of GUI.
%       Comments displayed at the command line in response 
%       to the help command. 
%

%
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

'test change'
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
          
%%  Construct the components
h=vPALM2_gui_construction(h,FIGH,FIGW);

%%  Initialization tasks

set(h.MainFig,'Visible','on')

set(h.MainFig, 'UserData', h);

%% Callbacks for MYGUI

h.h_vis=vPALM2_2d3d_visualization(h);

h.h_sliders=vPALM2_sliders(h);

set(h.MainFig,'CloseRequestFcn',{@vPALM2_newclosefunction,h});

%% ends
setcallbacks(h);

disp('GUI initialisation OK');


%%closing down
function vPALM2_newclosefunction(hObject, eventdata,h)
disp('everything is closing down!');

delete(h.h_vis.fig);
delete(h.h_sliders.fig);
delete(h.MainFigure3D);
delete(h.MainFig);

    

