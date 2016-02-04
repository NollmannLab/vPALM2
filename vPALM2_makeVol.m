function vPALM2_makeVol(x,y,z,VoxelSize,Filter,Treshold,handles)
% nBins(~150) defines the number of voxels in each direction. Be carefull not to use too high
%values ortherwise you'll get OUT OF MEMORY error messages.

% Filter=0 if you don't want 3D gaussian filtering, 1 if you do.

% Treshold=3 removes voxels with less than 3 localizations in it. Set it to
%0 if you don't want any thresholding
PixSize=str2num(get(handles.pixelsize,'String'));                                                                % 1pixel=100nm

list=get(handles.display3D,'String');
switch list{get(handles.display3D,'Value')}
    case 'scatter density'
        
            x=x-min(x);
            y=y-min(y);
            z=z-min(z);

            X=[x*PixSize,y*PixSize,z];

            radius=VoxelSize;

            D2=pdist2(X,X);

            D2(D2>radius)=0;
            D2(D2>0)=1;

            D2sum=sum(D2);

            figure('Color','k')
            scatter3(x*PixSize,y*PixSize,z,20,D2sum,'filled')
            set(gca,'Xcolor',[0.5 0.5 0.5]);
            set(gca,'Ycolor',[0.5 0.5 0.5]);
            set(gca,'Zcolor',[0.5 0.5 0.5]);
            set(gca,'Color',[0 0 0]);
            axis tight
            colormap hot
            colorbar 
        
    case 'voxels'

            % x and y values are in pixels units. z should be in 'nm' units.
            XZVoxelRatio=1;                                                           % Aspect ratio X/Z of the voxel
            nBinsX=round((max(x)-min(x))/VoxelSize*100);
            fprintf('Assuming a pixel size of %4.1f nm :\n',PixSize);
            fprintf('Voxels size along X axis are %4.1f nm\n', round((max(x)-min(x))/nBinsX*100));
            nBinsY= round ( nBinsX * round((max(y)-min(y)))/round((max(x)-min(x))));
            fprintf('Voxels size along Y axis are %4.1f nm\n', round((max(y)-min(y))/nBinsY*100));
            nBinsZ= round (nBinsX*  round((max(z)-min(z)))  / PixSize  /   round((max(x)-min(x)))  );
            fprintf('Voxels size along Z axis are %4.1f nm\n', round((max(z)-min(z))/nBinsZ*XZVoxelRatio));
            fprintf('The total number of voxels used is %8.0f out of ~6 000 000 (depending on the computer)\n',nBinsX*nBinsY*nBinsZ);

            xBins=linspace(min(x),max(x),nBinsX);                                       % Initialize Voxels arrays
            yBins=linspace(min(y),max(y),nBinsY);
            zBins=linspace(min(z),max(z),nBinsZ);
            D=zeros(nBinsX,nBinsY,nBinsZ);                                              % Initialize Voxels Matrix
            % D=zeros(nBins,nBins,nBins,'uint8'); 

            for i=1:numel(x)                                                            % Fill Voxels Matrix
            xi=find((x(i)> xBins),1,'last');
            yi=find((y(i)> yBins),1,'last');
            zi=find((z(i)> zBins),1,'last');
            D(xi,yi,zi)=D(xi,yi,zi)+1;
            end

            clear xBins yBins zBins xi yi zi

            if Filter>0
            h = fspecial3('gaussian',[2 2 4]);                                       % Filter parameters in X, Y, Z
            D=imfilter(D, h);                                                        % 3D filtering
            end
            D(D<Treshold)=0;                                                            % Treshold on the density voxels


            FIGW= 300;
            FIGH= 300;

            figure('Color','k');

            h = vPALM2_vol3d('cdata',D,'texture','3D');                                        % Plot the 3D density map


            %% Sets axis labels in nanometers
            NumTicks = 5;
            L = get(gca,'XLim');
            set(gca,'XTick',linspace(L(1),L(2),NumTicks))
            Xmin=round(min(x))*PixSize;
            Xmax=round(max(x))*PixSize;
            set (gca,'xTickLabel',0:((Xmax-Xmin)/(NumTicks-1)):(Xmax-Xmin))
            xlabel('X (nm)')

            L = get(gca,'YLim');
            set(gca,'YTick',linspace(L(1),L(2),NumTicks))
            Ymin=round(min(y))*PixSize;
            Ymax=round(max(y))*PixSize;
            set (gca,'yTickLabel',0:((Ymax-Ymin)/(NumTicks-1)):(Ymax-Ymin))
            ylabel('Y (nm)')

            L = get(gca,'ZLim');
            set(gca,'ZTick',linspace(L(1),L(2),NumTicks))
            Zmin=round(min(z));
            Zmax=round(max(z));
            set (gca,'zTickLabel',0:((Zmax-Zmin)/(NumTicks-1)):(Zmax-Zmin))

            zlabel('Z (nm)')

            axis tight;  daspect([1 ((Xmax-Xmin)/(Ymax-Ymin)) 0.25])
            grid on
            colorbar
            colormap hot
            set(gca,'Xcolor',[0.5 0.5 0.5]);
            set(gca,'Ycolor',[0.5 0.5 0.5]);
            set(gca,'Zcolor',[0.5 0.5 0.5]);
            set(gca,'Color',[0 0 0]);


            view(3)

    case 'surface'
        
          % x and y values are in pixels units. z should be in 'nm' units.
            XZVoxelRatio=1;                                                           % Aspect ratio X/Z of the voxel
            nBinsX=round((max(x)-min(x))/VoxelSize*100);
            fprintf('Assuming a pixel size of %4.1f nm :\n',PixSize);
            fprintf('Voxels size along X axis are %4.1f nm\n', round((max(x)-min(x))/nBinsX*100));
            nBinsY= round ( nBinsX * round((max(y)-min(y)))/round((max(x)-min(x))));
            fprintf('Voxels size along Y axis are %4.1f nm\n', round((max(y)-min(y))/nBinsY*100));
            nBinsZ= round (nBinsX*  round((max(z)-min(z)))  / PixSize  /   round((max(x)-min(x)))  );
            fprintf('Voxels size along Z axis are %4.1f nm\n', round((max(z)-min(z))/nBinsZ*XZVoxelRatio));
            fprintf('The total number of voxels used is %8.0f out of ~6 000 000 (depending on the computer)\n',nBinsX*nBinsY*nBinsZ);

            xBins=linspace(min(x),max(x),nBinsX);                                       % Initialize Voxels arrays
            yBins=linspace(min(y),max(y),nBinsY);
            zBins=linspace(min(z),max(z),nBinsZ);
            D=zeros(nBinsX,nBinsY,nBinsZ);                                              % Initialize Voxels Matrix
            % D=zeros(nBins,nBins,nBins,'uint8'); 

            for i=1:numel(x)                                                            % Fill Voxels Matrix
            xi=find((x(i)> xBins),1,'last');
            yi=find((y(i)> yBins),1,'last');
            zi=find((z(i)> zBins),1,'last');
            D(xi,yi,zi)=D(xi,yi,zi)+1;
            end

            clear xBins yBins zBins xi yi zi

            if Filter>0
            h = fspecial3('gaussian',[2 2 4]);                                       % Filter parameters in X, Y, Z
            D=imfilter(D, h);                                                        % 3D filtering
            end
            D(D<Treshold)=0;  
            
            figure, 
        
            cv = 5;
            Color = 'red';

            p=patch(isosurface(D,cv));
            set(p,'FaceColor',Color,'EdgeColor','none','FaceAlpha',.75);

            alpha(.8);
            daspect([1 1 1.])

            view(3)
            camlight left
            lighting gouraud

            
end



