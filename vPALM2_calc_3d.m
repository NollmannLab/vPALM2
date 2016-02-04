function m2=vPALM2_calc_3d(h)

wxPoly = h.wxPoly ;
wyPoly = h.wyPoly ;
ccPolyWtoZ = h.ccPolyWtoZ ;

nm_per_px = str2num(get(h.pixelsize,'String')); % nm
% step = 100 ; % nm

% automatically determine degree of polynomial fit
degree = 3% length(wxPoly) - 1 ;

% data


% [dataFilename,dataPathname] = uigetfile( '*.mat','Select MTT3D-generated .mat data file' ) ;
% load( [dataPathname dataFilename] ) ;

% .mat file organization
% row 1 : image (frame) number
% row 2 : x-coordinate
% row 3 : y-coordinate
% row 4 : alpha
% row 5 : -
% row 6 : 2D gaussian width
% row 7 : y-width
% row 8 : x-width
clear m;
m(:,1)=h.m(:,1); % frame number
m(:,2)=h.m(:,2); % pos x in px
m(:,3)=h.m(:,3); %pos y in px
m(:,4)=h.m(:,4); % intensity
m(:,5)=h.m(:,5); % resolution, nm
m(:,6)=h.m(:,13); %Wx
m(:,7)=h.m(:,13); %Wx
m(:,8)=h.m(:,14); %Wy

m(:,10)=h.m(:,10); %angle
m(:,11)=h.m(:,11); %asymmetry
m(:,12)=h.m(:,12); %channel
m(:,13)=h.m(:,8); %slice
m(:,9)=h.m(:,9); %position

matrice_results=m';
% remove zero columns ;
% matrice_results( :,matrice_results(1,:)==0 ) = [] ;

% remove overly high intensities
% matrice_results( :,matrice_results(4,:)>10^5 ) = [] ;

% plotting data matrix
% column 1 : x-coordinate
% column 2 : y-coordinate
% column 3 : z-coordinate
% column 4 : x-width
% column 5 : y-width
% column 6 : x-width minus y-width
% column 7 : intensity
% column 8 : frame

plotData = zeros( length(matrice_results),8 ) ;
plotData( :,1 ) = [matrice_results( 2,: )*nm_per_px]' ; % x-coordinate
plotData( :,2 ) = [matrice_results( 3,: )*nm_per_px]' ; % y-coordinate
plotData( :,4 ) = matrice_results( 7,: )' ; % y-width 
plotData( :,5 ) = matrice_results( 8,: )' ; % x-width
plotData( :,6 ) = matrice_results( 8,: )' - matrice_results( 7,: )' ; % wx-wy
plotData( :,7 ) = matrice_results( 4,: )' ; % intensity
plotData( :,8 ) = matrice_results( 1,: )' ; % frame
plotData( :,9 ) = matrice_results( 5,: )' ; % resolution
plotData( :,10 ) = matrice_results( 6,: )' ; % bck
plotData( :,11 ) = matrice_results( 10,: )' ; % angle
plotData( :,12 ) = matrice_results( 11,: )' ; % asymmetry
plotData( :,13 ) = matrice_results( 12,: )' ; % channel
plotData( :,14 ) = matrice_results( 13,: )' ; %slice
plotData( :,15 ) = matrice_results( 9,: )' ; %position

% determine z-coordinates based on wx-wy
curveThreshold = 1 ; % px
removedPoints = 0 ;
for q = 1:1:length(matrice_results)
	plotData(q,3) = polyval( ccPolyWtoZ,plotData(q,6) ) ;
    % filter out points curveThreshold [nm] away from wx and wy curves
	idealWidthX = polyval(wxPoly,plotData(q,3)) ;
	idealWidthY = polyval(wyPoly,plotData(q,3)) ;
	if ( (abs(idealWidthY-plotData(q,4)) > curveThreshold) || (abs(idealWidthX-plotData(q,5)) > curveThreshold) )
		plotData(q,3) = nan ;
		removedPoints = removedPoints + 1 ;
	end
end

% remove filtered out data
plotData( isnan(plotData(:,3)),: ) = [] ;
% plotData( plotData(:,3)<-1000,: ) = [] ;
% plotData( plotData(:,3)>1000,: ) = [] ;

% find max and min of delta calibration curve
[maxZ,maxZIndex] = max( plotData(:,3) ) ;
[minZ,minZIndex] = min( plotData(:,3) ) ;

% plotData( plotData(:,6) > plotData(maxZIndex,6) | plotData(:,6) < plotData(minZIndex,6),: ) = []; 

D = [ zeros(length(plotData(:,1)),1) plotData(:,1) plotData(:,2) plotData(:,3) plotData(:,8) plotData(:,7) ] ;
D = sortrows( D,-6 ) ;

m2(:,1)=plotData(:,8); %fr
m2(:,2)=plotData(:,1)./nm_per_px; %x
m2(:,3)=plotData(:,2)./nm_per_px; %y
m2(:,4)=plotData(:,7); %I
m2(:,5)=plotData(:,9); % resolution
m2(:,6)=plotData(:,10); % bck
m2(:,8)=plotData(:,14);% slice
m2(:,9)=plotData(:,15);% position
m2(:,10)=plotData(:,11); % angle
m2(:,11)=plotData(:,12); % asymmetry
m2(:,12)=plotData(:,13);% channel
m2(:,13)=plotData(:,4); %wx
m2(:,14)=plotData(:,5); %wy
m2(:,15)=plotData(:,3); %z


% add index column
for u = 1:1:(length(plotData))
	D(u,1) = u ;
end		

deltaRange = min(plotData(:,6)):0.02:max(plotData(:,6)) ;
zRange = min(plotData(:,3)):0.02:max(plotData(:,3)) ;

% plot raw data on wx-wy curve
% plot to show width distribution
figure,
subplot(2,1,1) ;
scatter( plotData(:,3),plotData(:,6),'o','markeredgecolor',[0 0 1] ) ; hold on ;
% scatter( plotData(:,8)*step-alignmentOffset,plotData(:,6),'.','markeredgecolor',[0 0 1] ) ; 
plot( polyval(ccPolyWtoZ,deltaRange),deltaRange,'linewidth',1,'color','k' ) ; hold on ;
ylabel( 'w_x-w_y [px]','fontsize',16 ) ;
xlabel( 'z-Position [nm]','fontsize',16 ) ;
title( h.fullFileName,'fontsize',16,'Interpreter','none' ) ;
set(get(gcf,'CurrentAxes'),'FontSize',14) ;
grid ;
% xlim([-1000 1000])
legend( 'Raw Data Points','Calibration Curve','BestOutside' ) ;

% plot to show width distribution

subplot(2,1,2) ;
scatter( plotData(:,3),plotData(:,4),'.','markeredgecolor',[0 0 1] ) ; hold on ;
scatter( plotData(:,3),plotData(:,5),'.','markeredgecolor',[0 1 0] ) ; hold on ;
plot( zRange,polyval(wxPoly,zRange),'linewidth',1,'color','k' ) ; hold on ; 
plot( zRange,polyval(wyPoly,zRange),'linewidth',1,'color','k' ) ; hold on ;
ylabel( 'PSF Width [px]','fontsize',16 ) ;
xlabel( 'z-Position [nm]','fontsize',16 ) ;
title( h.fullFileName,'fontsize',16,'Interpreter','none' ) ;
% ylim( [0 3.5] ) ;
% xlim([-1000 1000])
set(get(gcf,'CurrentAxes'),'FontSize',14) ;
grid ;
legend( 'w_x Fit Data','w_y Fit Data','w_x Curve','w_y Curve','BestOutside' ,'fontsize',10) ;

% [filename,pathname] = uiputfile( '*.roi','Save ROI',dataFilename ) ;

% header string
% header = [ 'index' 'x [nm]' 'y [nm]' 'z [nm]' 'frame' 'intensity' ] ;
% dlmwrite( [pathname filename],header,'delimiter','\t' ) ;
% dlmwrite( [pathname filename],D,'delimiter','\t','-append' ) ;

