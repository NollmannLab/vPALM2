function [wxPoly,wyPoly,ccPolyWtoZ]=vPALM2_cal3D(h)


step = str2num(get(h.stepsize,'String'));%i//nput( 'Enter step size of z-stack [nm]: ' ) ;
degree = 3;%input( 'Select degree of fit [3-7]: ') ; 

%% .mat file organization
%  row 1 : image (frame) number
%  row 2 : x-coordinate
%  row 3 : y-coordinate
%  row 4 : alpha
%  row 5 : -
%  row 6 : 2D gaussian width
%  row 7 : x-width
%  row 8 : y-width

if ~ischar(h.fullFileName)
	files = length(h.fullFileName) ;
else
	files = 1 ;
end

polyVecsX = zeros( files,degree+1 ) ;
polyVecsY = zeros( files,degree+1 ) ;
polyVecsDelta = zeros( files,degree+1 ) ;
zPositions = zeros( files,600 ) ; % assume frames vector length < 600

% store cutoff positions
cutoffsX = zeros( files,2 ) ;
cutoffsY = zeros( files,2 ) ;

frameData = [] ;
zPositionData = [] ;
wyData = [] ;
wxData = [] ;
deltaData = [] ;

figure
	for k = 1:1:files

m(:,1)=h.m(:,1); % frame number
m(:,2)=h.m(:,2); % pos x in px
m(:,3)=h.m(:,3); %pos y in px
m(:,4)=h.m(:,4); % intensity
m(:,5)=h.m(:,5); % resolution, nm
m(:,6)=h.m(:,13); %Wx
m(:,7)=h.m(:,13); %Wx
m(:,8)=h.m(:,14); %Wy
matrice_results=m';

	% remove zero columns ;
	matrice_results( :,matrice_results(1,:)==0 ) = [] ;
	% remove columns with wx or wy greater than 3.5 px
	matrice_results( :,matrice_results(7,:)>10 ) = [] ; %3
	matrice_results( :,matrice_results(8,:)>10 ) = [] ; %3
	matrice_results( :,matrice_results(7,:)<0.5 ) = [] ;
	matrice_results( :,matrice_results(8,:)<0.5 ) = [] ;	

	% select frame cutoff values
	subplot(2,1,1)
	plot( matrice_results(1,:),matrice_results(8,:),'o-b' ) ; hold on ;
	plot( matrice_results(1,:),matrice_results(7,:),'o-r' ) ; 
	grid on ;
	xlabel( 'Frame Number','fontsize',16 ) ;
	ylabel( 'x,y PSF Widths [px]','fontsize',16 ) ;
	set(get(gcf,'CurrentAxes'),'FontSize',14) ;
	legend('w_x','w_y','location','Best') ;
	ylim([0 10]) ;
% 	if (files>1)
% 		title( current_filename ) ;
% 	else
% 		title( filename ) ;
% 	end
	k=1;
	subplot(2,1,2)
	plot( matrice_results(1,:),matrice_results(8,:)-matrice_results(7,:),'o-b' ) ; hold on ;
	grid on ;
	xlabel( 'Frame Number','fontsize',16 ) ;
	ylabel( 'w_x - w_y [px]','fontsize',16 ) ;
	set(get(gcf,'CurrentAxes'),'FontSize',14) ;
	ylim([-2.5 2.5]) ;
    fr=matrice_results(1,:);
% 	lowerCutoffX = input( 'Enter lower frame cutoff: ' ) ;
% 	upperCutoffX = input( 'Enter upper frame cutoff: ' ) ;
    disp('Click first point of event to analyze:')
    regionstart= ginput(1);
    regionstart= find(fr> regionstart(1) );
    regionstart= min(regionstart);
    lowerCutoffX= fr(regionstart);
    disp('Click last point of event to analyze:')
    regionend= ginput(1);
    regionend= find(fr < regionend(1) );
    regionend= max(regionend);
    upperCutoffX= fr(regionend);   
    
	lowerCutoffY = lowerCutoffX ;
	upperCutoffY = upperCutoffX ;
	
% 	lowerCutoffY = input( 'Enter lower frame y cutoff: ' ) ;
% 	upperCutoffY = input( 'Enter upper frame y cutoff: ' ) ;

	cutoffsX( k,1 ) = lowerCutoffX*step ;
	cutoffsX( k,2 ) = upperCutoffX*step ;

	cutoffsY( k,1 ) = lowerCutoffY*step ;
	cutoffsY( k,2 ) = upperCutoffY*step ;
	
	% x-points removal
	matrice_results( 8,matrice_results(1,:)<lowerCutoffX ) = nan ;
	matrice_results( 8,matrice_results(1,:)>upperCutoffX ) = nan ;

	% y-points removal
	matrice_results( 7,matrice_results(1,:)<lowerCutoffY ) = nan ;
	matrice_results( 7,matrice_results(1,:)>upperCutoffY ) = nan ;

	% frame removal
	matrice_results( 1,matrice_results(1,:)<lowerCutoffY ) = nan ;
	matrice_results( 1,matrice_results(1,:)>upperCutoffY ) = nan ;
	
	frame = matrice_results( 1,: )' ;
	zPosition = frame*step ;
	wy = matrice_results( 7,: )' ;
	wx = matrice_results( 8,: )' ;
	delta = matrice_results( 8,: )' - matrice_results( 7,: )' ;

	% store data from each bead in column vector
	if (k > 1)
		dataDiff = length(frameData) - length(frame) ;
		if (dataDiff < 0)
			frameData = [ frameData ; nan*zeros(abs(dataDiff),k-1) ] ;
			zPositionData = [ zPositionData ; nan*zeros(abs(dataDiff),k-1) ] ;
			wyData = [ wyData ; nan*zeros(abs(dataDiff),k-1) ] ;
			wxData = [ wxData ; nan*zeros(abs(dataDiff),k-1) ] ;
			deltaData = [ deltaData ; nan*zeros(abs(dataDiff),k-1) ] ;
		elseif (dataDiff > 0)
			frame = [ frame ; nan*zeros(dataDiff,1) ] ;
			zPosition = [ zPosition ; nan*zeros(dataDiff,1) ] ;
			wy = [ wy ; nan*zeros(dataDiff,1) ] ;
			wx = [ wx ; nan*zeros(dataDiff,1) ] ;
			delta = [ delta ; nan*zeros(dataDiff,1) ] ;
		end
	end

	frameData = [frameData , frame] ;
	zPositionData = [zPositionData , zPosition] ;
	wyData = [wyData , wy] ;
	wxData = [wxData , wx] ;
	deltaData = [deltaData , delta] ;
	
	% find beginning and end of wx,wy vectors
	beginX = nan ;
	t = 1 ;
	while (isnan(beginX))
		if ~isnan(wx(t)) ;
			beginX = t ;
		end
		t = t + 1 ;
	end
	
	endX = nan ;
	t = length(frame) ;
	while (isnan(endX))
		if ~isnan(wx(t)) ;
			endX = t ;
		end
		t = t - 1 ;
	end
	
	beginY = nan ;
	t = 1 ;
	while (isnan(beginY))
		if ~isnan(wy(t)) ;
			beginY = t ;
		end
		t = t + 1 ;
	end
	
	endY = nan ;
	t = length(frame) ;
	while (isnan(endY))
		if ~isnan(wy(t)) ;
			endY = t ;
		end
		t = t - 1 ;
	end

	% fit
	polyVecsX( k,: ) = polyfit( zPosition(beginX:endX),wx(beginX:endX),degree ) ;
	polyVecsY( k,: ) = polyfit( zPosition(beginY:endY),wy(beginY:endY),degree ) ;
	polyVecsDelta( k,: ) = polyfit( zPosition(beginY:endX),delta(beginY:endX),degree ) ;

	polyVecsXplot( 1,: ) = polyfit( frame(beginX:endX),wx(beginX:endX),degree ) ;
	polyVecsYplot( 1,: ) = polyfit( frame(beginY:endY),wy(beginY:endY),degree ) ;
	polyVecsDeltaPlot( 1,: ) = polyfit( frame(beginY:endX),delta(beginY:endX),degree ) ;

	% store z-positions
	zPositions( k,1:length(zPosition) ) = zPosition ;

	% overlay fit to raw data
	polyFitXplot = polyval( polyVecsXplot,frame ) ;
	polyFitYplot = polyval( polyVecsYplot,frame ) ;
	polyFitDeltaPlot = polyval( polyVecsDeltaPlot,frame ) ;

	subplot(2,1,1)
	plot( frame,polyFitXplot,'k' ,'Linewidth',2) ; hold on ;
	plot( frame,polyFitYplot,'k' ,'Linewidth',2) ; hold off ;
	
	subplot(2,1,2)
	plot( frame,polyFitDeltaPlot,'k','Linewidth',2 ) ;
	hold off ;
	
% 	disp('Press any button to continue ...')
% 	w = waitforbuttonpress;

% 	while (w~=1) ;
% 	end ;	
	

% store roots of polynomials for each file
polyRoots = zeros(files,1) ;

for b = 1:1:files

	diffPoly = polyVecsDelta( b,: ) ; % polyVecsX( b,: ) - polyVecsY( b,: ) ;
	diffRoots = roots(diffPoly) ;
	diffRoots = diffRoots(find(imag(diffRoots)==0)) ; % select only real roots

	for m = 1:1:length(diffRoots)
		if ( diffRoots(m) > cutoffsY( b,1) ) && ( diffRoots(m) < cutoffsX( b,2 ) ) 
			root = diffRoots(m) ;
		end
	end
end
    
	polyRoots( b ) = root ;
	
end

% adjust raw data with offset (root), i.e., centre raw data to z = 0
for c = 1:1:files ;
	zPositionData(:,c) = zPositionData(:,c) - polyRoots(c) ;
end

% reshape raw data vectors for plotting
zPositionData = reshape( zPositionData,numel(zPositionData),1 )  ;
deltaData = reshape( deltaData,numel(deltaData),1 )  ;
wxData = reshape( wxData,numel(wxData),1 )  ;
wyData = reshape( wyData,numel(wyData),1 )  ;

% remove NaN values
zPositionData(isnan(zPositionData(:))) = [] ;
deltaData(isnan(deltaData(:))) = [] ;
wxData(isnan(wxData(:))) = [] ;
wyData(isnan(wyData(:))) = [] ;

% sort data
allData = [zPositionData deltaData wxData wyData] ;
allData = sortrows(allData,1) ; % sort based on zPosition
zPositionData = allData( :,1 ) ;
deltaData = allData( :,2 ) ;
wxData = allData( :,3 ) ;
wyData = allData( :,4 ) ;

% fit
wxFit = polyfit( zPositionData,wxData,degree ) ;
wyFit = polyfit( zPositionData,wyData,degree ) ;
deltaFit = polyfit( zPositionData,deltaData,degree ) ;

% overlay fit to raw data
wxFitPlot = polyval( wxFit,zPositionData ) ;
wyFitPlot = polyval( wyFit,zPositionData ) ;
deltaFitPlot = polyval( deltaFit,zPositionData ) ;

figure
for k = 1:1:files ;
	subplot(2,1,1) ;
	scatter( zPositionData,deltaData,'o','filled','markeredgecolor',[0 0 0],'MarkerFaceColor',[0 0 1] ) ; hold on ;
	plot( zPositionData,deltaFitPlot,'k','linewidth',1.5 ) ; hold on ;
	subplot(2,1,2) ;
    scatter( zPositionData,wxData,'o','filled','markeredgecolor',[0 0 0],'MarkerFaceColor',[0 0 1] ) ; hold on ;
    scatter( zPositionData,wyData,'o','filled','markeredgecolor',[0 0 0],'MarkerFaceColor',[1 0 0] ) ; hold on ;
    plot( zPositionData,wxFitPlot,'k','linewidth',1.5 ) ; hold on ;
    plot( zPositionData,wyFitPlot,'k','linewidth',1.5 ) ; hold on ;
end	

subplot(2,1,1) ;
grid on ;
% xlabel( 'z-Position [nm]','fontsize',16 ) ;
ylabel( 'w_x-w_y [px]','fontsize',16 ) ;
set(get(gcf,'CurrentAxes'),'FontSize',14) ;
% legend('w_x-w_y','Fit','location','Best') ;
ylim([-2 2]) ;
xlim([-400 400])
legend1=legend('w_x-w_y','Fit','location','Best');
set(legend1,'FontSize',11,'Box','off');

subplot(2,1,2) ;
grid on ;
xlabel( 'z-Position [nm]','fontsize',16 ) ;
ylabel( 'w_x, w_y [px]','fontsize',16 ) ;
set(get(gcf,'CurrentAxes'),'FontSize',14) ;
% legend('w_y Data','w_x Data','Fit','location','Best') ;
ylim([0.5 10]) ;
xlim([-600 600])
legend2=legend('w_y Data','w_x Data','Fit','location','Best');
set(legend2,'FontSize',11,'Box','off');

% condition polynomials for output
wxPoly = polyfit(zPositionData,wxFitPlot,degree) ;
wyPoly = polyfit(zPositionData,wyFitPlot,degree) ;
ccPolyWtoZ = polyfit(deltaFitPlot,zPositionData,degree) ;

%% output calibration curve
% format
% row 1 : averaged wx polynomial coefficients
% row 2 : averaged wy polynomial coefficients
% row 3 : averaged wx-wy polynomial coefficients
% row 4 : alignment offset

% [fileTitle,pathnameTitle] = uiputfile( '*.clb', 'Save calibration file' ) ;
% 
% dlmwrite( [pathnameTitle fileTitle],['wx polynomial:	' num2str(wxFit,'%0.16e\t')],'delimiter','' ) ;
% dlmwrite( [pathnameTitle fileTitle],['wy polynomial:	' num2str(wyFit,'%0.16e\t')],'delimiter','','-append' ) ;
% dlmwrite( [pathnameTitle fileTitle],['dw polynomial:	' num2str(ccPolyWtoZ,'%0.16e\t')],'delimiter','','-append' ) ;
