function h = vPALM2_Beads_detection(h)

q =  str2num(get(h.pixelsize,'String'));

%% 1- Load an image for the beads selection (easier to use an image than a
%% plot of all the localizations) and detect the brightest beads.

[ImageName,ImagePath] = uigetfile('*.tif','Open image for crop selection');
FullImageName = strcat(ImagePath,ImageName);
tiff_info=imfinfo(FullImageName); %Return a structure containing all the information for EACH image of the movie
Nframes=length(tiff_info); %Return the number of frames
im = imread(FullImageName,round(Nframes/2)); % Display the image
im = imadjust(im,[0 0.1],[]); % Adjust the image contrast

fig1 = figure;
set(0,'Units','pixels'); %Define the type of units used later for the position (here in pixels)
scnsize = get(0,'ScreenSize');%Get the size of the screen in pixels
set(fig1,'OuterPosition',scnsize);%Display fig1 in order to completely fill the screen

hfig1 = imagesc(im);
hold on
colormap('gray');
axis image

Max = max(max(im));
Min = min(min(im));
tiff_info = tiff_info(1);
Width = tiff_info.Width;
Height = tiff_info.Height;

Binary = zeros(Height,Width);
Level = round(Min + (Max-Min)*0.5);

for k = 1 : Width
    for n = 1 : Height
        if im(n,k)>Level
            Binary(n,k) = 1;
        end
    end
end

Binary = imfill(Binary,'holes');
Boundaries = bwboundaries(Binary);
NBoundaries = length(Boundaries);
ROIs = {};

for k = 1:NBoundaries
    boundary = Boundaries{k};
    X = round(mean(boundary(:,2)));
    Y = round(mean(boundary(:,1)));
    rectangle('Position',[X-10,Y-10,20,20],'EdgeColor','g')
    ROIs{k} = [X-10,Y-10,20,20];
end

%% 2- Allow for the manual detection of beads and add all the ROIs to
%% the structure ROIs

NROIsTotal = NBoundaries;
Manual_Detection = questdlg('Do you want to select other beads?','Manual Selection','Yes','No','Cancel','Yes');

switch Manual_Detection
    
    case 'Yes'
        
        hSP = imscrollpanel(fig1,hfig1);
        set(hSP,'HitTest','off');
        api = iptgetapi(hSP);
        api.setMagnification(4) % 2X = 200%
        fig_overview = imoverview(hfig1);
        
        Proceed = 1;
        while Proceed
            
            NROIsTotal = NROIsTotal + 1;
            rect = getrect(fig1);
            
            Accept_rectangle = 0;
            ProceedTo_Next_ROI = 0;
            
            hrect(1) = imrect(gca,rect);
            AcceptROI = 0;
            
            while ~ProceedTo_Next_ROI
                
                k = waitforbuttonpress;
                if k % Keyboard press
                    KeyPress = get(gcf,'CurrentCharacter');
                    switch KeyPress
                        case ' '
                            ProceedTo_Next_ROI = 1;
                            AcceptROI = 1;
                        case 'f'
                            ProceedTo_Next_ROI = 1;
                            Proceed = 0;
                            AcceptROI = 1;
                    end
                else % Mouse press
                    st = get(gcf,'SelectionType');
                    switch st
                        case {'open'} % double click to remove a rectangle
                            if ~isempty(gco)
                                delete(hrect(1));
                                move2next_rect = 1;
                                NROIsTotal = NROIsTotal - 1;
                                ProceedTo_Next_ROI = 1;
                            end
                    end
                end
            end
            
            if AcceptROI
                api = iptgetapi(hrect(1));
                rect = api.getPosition();
                ROIs{NROIsTotal} = rect;
                rectangle('Position',rect,'EdgeColor','y','LineWidth',2,'LineStyle','--');
                text(rect(1)+5,rect(2),num2str(NROIsTotal),'Color','y','FontSize',12,'VerticalAlignment','Top');
                delete(hrect(1));
            end
            %             if ~Proceed
            %                 close(fig_overview)
            %
            %             end
            
        end
end

close(fig1)
axes(h.axes1), hold off
imagesc(im); hold on, colormap('gray'), axis image
for k = 1 : NROIsTotal
    rectangle('Position',ROIs{k},'EdgeColor','g');
    text(ROIs{k}(1)+5,ROIs{k}(2),num2str(k),'Color','y','FontSize',12,'VerticalAlignment','Top');
end


%% 3- Detection of the beads and calculation of the trajectories

M = h.m; % Load the events
M = M(:,1:3);

[~,Idx] = sort(M(:,1));
M = M(Idx,:);

switch get(h.load_options,'Value')
    case 4 %MTT
        M1(:,1)=M(:,1);
        M1(:,3)=M(:,2)+1;
        M1(:,2)=M(:,3)+1;
        M=M1;
        
    case 6 %micromanager
        
        ;
    case 7 %micromanager 3D
        ;
        
end

figure
imshow(im)
hold on
plot(M(:,2),M(:,3),'.b')

% figure
% imshow(im)
% hold on
% plot(M(:,3)+1,M(:,2)+1,'.b')

AnalysisMethods = get(h.load_options,'Value'); % Return the method used to detect the single particle fluorescent events
MaxFrame = M(end,1); % Return the maximum number of frames
Beads = struct;
% M(:,2)=M(:,2)+2;
% M(:,3)=M(:,3)+6;

hold on
Test = zeros(NROIsTotal,MaxFrame);

for nROI = 1 : NROIsTotal
    Position{nROI} = zeros(1,MaxFrame);
end
%
% switch AnalysisMethods
%
%     case 3

hwb = waitbar(0,'Calculating the trajectories...');
MaxLine = length(M);

for Frame = 1 : MaxFrame
    
    if Frame == 1
        Line = 1;
        LineStart = 1;
        while M(Line,1) == 1
            Line = Line +1;
        end
        LineStop = Line - 1;
    else
        LineStart = LineStop + 1;
        while M(Line,1) == Frame
            Line = Line +1;
            if Line>MaxLine
                break
            end
        end
        LineStop = Line - 1;
    end

    number_frame = LineStart:1:LineStop;
    %             m_frame = M(number_frame,:); % Look for all the events that have been detected on the image # 'Frame'
    %             m_frame = cat(2,M(number_frame,1),M(number_frame,3),512-M(number_frame,2));
    m_frame = cat(2,M(number_frame,1),M(number_frame,2),M(number_frame,3));
    
    for nROI = 1 : NROIsTotal
        
        waitbar(Frame/MaxFrame); % Give an indication of how fast the process is going
        ROI = ROIs{nROI};
        Xmin = round(ROI(1,1));
        Ymin = round(ROI(1,2));
        Xmax = round(ROI(1,1)+ROI(1,3));
        Ymax = round(ROI(1,2)+ROI(1,4));
        
        Events_ROI = find(m_frame(:,2)<=Xmax & m_frame(:,2)>=Xmin & m_frame(:,3)<Ymax & m_frame(:,3)>=Ymin);
        %                 Events_ROI = find(m_frame(:,2)<=Ymax & m_frame(:,2)>=Ymin & m_frame(:,3)<Xmax & m_frame(:,3)>=Xmin);
        
        Test(nROI,Frame) = size(Events_ROI,1);
        
        if size(Events_ROI,1) == 1
            Position{nROI}(1,Frame) = number_frame(Events_ROI); % Return the position in M of the bead detected in the ROI
        else
%             Frame
%             size(Events_ROI,1)
%             Position{nROI}(:,Frame) = [];
        end
    end
end

for nROI = 1 : NROIsTotal
    Idx = find(Position{nROI}(1,:)==0);
    Position{nROI}(:,Idx) = [];
end

for nROI = 1 : NROIsTotal
    NameField = strcat('Bead_',num2str(nROI));
    LBead(nROI) = length(Position{nROI});
    Beads.(NameField) = M(Position{nROI},:);
end


%% 4- Analysis of the trajectories. The trajectories that have been validated are analysed here and saved
%% in the structures called "trajectories" and "trajectories-sg".

NRejected = 0;

for Number = 1 : NROIsTotal
    
    waitbar((Number-NRejected)/NROIsTotal);
    clear bead_sg; %It needs to be done, else there is an error when the matrix do not have exactly the same size (but not sure I understand why...)
    clear bead;
    
    NameField = strcat('Bead_',num2str(Number-NRejected));
    bead = Beads.(NameField);
    
    if length(bead)<0.5*MaxFrame % If too few events are detected within the cropped area the bead is discarted
        warndlg(strcat('The bead #', num2str(Number), ' is discarted since there were too few events detected'),'!! Warning !!');
        uiwait;
        
        % The beads that have a number higher than "Number"
        % are renamed below.
        % -----------------
        LBead(Number-NRejected) = [];
        [NROIsTotal, Beads] = Rename_Beads(Number-NRejected, NROIsTotal, Beads);
        NRejected = NRejected + 1;
        
    else
        [bead, Accepted] = Bead_check_qPALM(bead', MaxFrame, Number-NRejected); % This function makes sure that for each frame there is one and only one event detected
        
        if isequal(Accepted, 'Yes')
            bead(1,:) = bead(1,:)*50/60000; % Transcribe the # of frames in time
            
            bead_sg(2,:) = smooth(bead(2,:),21,'sgolay',1); % smoothing
            bead_sg(3,:) = smooth(bead(3,:),21,'sgolay',1);
            
            if ~exist('bead0','var')
                bead0 = bead;
                ZeroX_all = mean(bead(2,1:15)) ;
                ZeroY_all = mean(bead(3,1:15)) ;
            else
                % The first bead accepetd is used a a reference and all the other
                % beads are shifted in such way that the distance between
                % both trajectories is the smallest possible.
                % ------------------------------------------
                
                dx = sum(bead0(2,:)-bead(2,:))/length(bead0(2,:));
                dy = sum(bead0(3,:)-bead(3,:))/length(bead0(3,:));
                bead_sg(2,:) = bead_sg(2,:) + dx;
                bead_sg(3,:) = bead_sg(3,:) + dy;
                ZeroX_all = [ZeroX_all , mean(bead(2,1:15)+dx)] ;
                ZeroY_all = [ZeroY_all , mean(bead(3,1:15)+dy)] ;
            end
            
        elseif ~isequal(Accepted, 'Yes')
            
            LBead(Number-NRejected) = [];
            [NROIsTotal, Beads] = Rename_Beads(Number-NRejected, NROIsTotal, Beads);
            NRejected = NRejected + 1;
            
        end
        % If the calculation went well and the bead have been accepted,
        % the trajectories are calculate here. Else, this part is simply
        % skipped
        % ------
        if exist('bead','var') && exist('bead_sg','var')
            frames = bead(1,:);
            X = bead(2,:);
            Y = bead(3,:);
            X_sg = bead_sg(2,:);
            Y_sg = bead_sg(3,:);
            trajectories.(NameField) = {frames, X, Y};
            trajectories_sg.(NameField) = {frames, X_sg, Y_sg};
        end
    end
end

close(hwb);

%% 5- Calculation of the origin.

% The lists DX and DY contain the means over the 20 first frames of EACH
% trajectories. mean(DX) and mean(DY) will be used to make sure that all
% trajectories are actually starting at zero, despite the noise.
% -----------------------------------------------------------

if exist('ZeroX_all','var') && exist('ZeroY_all','var')
    ZeroX = mean(ZeroX_all);
    ZeroY = mean(ZeroY_all);
    ZeroX_all = ZeroX_all - ZeroX;
    ZeroY_all = ZeroY_all - ZeroY;
else
    warndlg('There is no bead to correct the drift with. The calculation is stopped here.')
    uiwait
    return
end

for n = 1 : NROIsTotal
    
    waitbar(n/NROIsTotal)
    clear('bead', 'bead_sg', 'X', 'Y', 'X_sg', 'Y_sg');
    
    NameField = strcat('Bead_',num2str(n));
    bead_sg = trajectories_sg.(NameField);
    
    frames = bead_sg{1};
    X_sg = bead_sg{2};
    Y_sg = bead_sg{3};
    
    X_sg = X_sg - ZeroX;
    Y_sg = Y_sg - ZeroY;
    
    trajectories_sg.(NameField) =  {frames, X_sg, Y_sg};
end

PositionX = zeros(NROIsTotal);
PositionY = zeros(NROIsTotal);

%% 6- Display the trajectories. Allow the remove the trajectories that are clearly not good. Calculate the
%% reference trajectory.

Colors = [[1 0 0];[0 1 0];[0 0 0];[1 1 0];[1 0 1];[0 1 1];[0.5 0.5 0];...
    [0.5 0 0.5];[0 0.5 0.5];[0.5 0.5 0.5];[1 0.5 0];[0.5 1 0];[1 0 0.5];...
    [0.5 0 1];[0 1 0.5];[0 0.5 1];[1 0.5 0.5];[0.5 1 0.5];[0.5 0.5 1];...
    [0 0 0.25];[0.25 0 0];[0 0.25 0];[0 0.5 0.25];[0.25 0.25 0];...
    [0.25 0 0.25];[0 0.25 0.25]]; %Matrix that will help to draw the curves in different colours

if NROIsTotal>0
    
    Proceed = 0;
    clear Answer;
    while Proceed == 0
        
        %% Define the size of the figure with respect to the size of the screen
        
        fig1 = figure;
        set(0,'Units','pixels'); %Define the type of units used later for the position (here in pixels)
        scnsize = get(0,'ScreenSize');%Get the size of the screen in pixels
        set(fig1,'OuterPosition',scnsize);%Display fig1 in order to completely fill the screen
        
        %% For each bead, the X and Y trajectories are plotted as well as their positions
        %% on the image. A specific color is assigned for each bead, in order to make
        %% the analysis a little bit easier
        
        clear('PositionX', 'PositionY')
        for n = 1 : NROIsTotal
            
            if n > size(Colors,1)
                k = n - floor(n/size(Colors,1))*size(Colors,1);
            else
                k = n;
            end
            
            NameField = strcat('Bead_',num2str(n));
            
            T = trajectories.(NameField);
            T_sg = trajectories_sg.(NameField);
            Frames = T_sg{1,1};
            X = T{1,2};
            Y = T{1,3};
            X_sg = T_sg{1,2};
            Y_sg = T_sg{1,3};
            
            subplot(1,3,1)
            hold on
            
            plot(Frames,X_sg,'Color',Colors(k,:))
            xlabel('time(min)')
            ylabel('drift along x')
            
            subplot(1,3,2)
            hold on
            
            plot(Frames,Y_sg,'Color',Colors(k,:))
            xlabel('time(min)')
            ylabel('drift along y')
            
            subplot(1,3,3)
            hold on
            
            if n==1
                imagesc(im); colormap('gray'), axis image
            end
            
            X = X(1:10)+ 1;
            Y = Y(1:10)+ 1;
            Xmax = max(X);
            Ymax = max(Y);
            PositionX(n) = mean(X); % Arrays that will keep in memory the X and Y positions of the selected beads on the image
            PositionY(n) = mean(Y);
            
            plot(X,Y,'--rs','LineWidth',2, 'Color',Colors(k,:),'MarkerSize',10)
            text(Xmax+7,Ymax+7,num2str(n),'FontSize',10,'Color',Colors(k,:))
            title('Position of each selected bead on the image')
            axis square
        end
        
        subplot(1,3,1)
        axis tight
        subplot(1,3,2)
        axis tight
        
        %% Quite often, there one or more beads that are moving too much to allow for a proper
        %% correction of the drift. They can be removed manually by the user. The question is asked only
        %% once. If the answer is 'Yes', the function will run as long as there is at least one bead left
        %% or the user is not cliking right.
        
        if ~exist('Answer','var')
            Answer  = questdlg('Is there a bead to remove? (POINT IT ON THE IMAGE, NOT ON THE GRAPHS)','Correction','Yes','No','Yes');
        end
        switch Answer
            case 'Yes'
                
                Xsubplot = [0.13 ; 0.41 ; 0.69]; %Position of the graphs along the X axis
                
                keepgoing = 0;
                while ~keepgoing
                    
                    Button_press = waitforbuttonpress;
                    if Button_press % Means the user pressed a key on the keyboard
                        Stop_Removing_Bead = questdlg('Are you done removing beads?','Done?','Yes','No','Yes');
                        switch Stop_Removing_Bead
                            case 'Yes'
                                Proceed = 1;
                                break
                        end
                    elseif ~Button_press % Means the mouse was used
                        
                        Right_Left = get(gcf,'SelectionType'); %Returns if it was Right or Left
                        switch Right_Left
                            
                            case 'normal'
                                
                                C_gcf = get (gcf, 'CurrentPoint'); %Return the position of the mouse when the user click on it (for the while figure)
                                C_gca = get (gca, 'CurrentPoint'); %Return the position of the mouse when the user click on it (for the current axis system)
                                X = C_gcf(1,1); %Calculate the relative X position of the mouse
                                [~,nX] = min(abs(Xsubplot-X/scnsize(1,3))); % Look for the X position of the closest subplot
                                N = nX;
                                
                                if N == 3
                                    D = sqrt((PositionX(:)-C_gca(1,1)).^2 + (PositionY(:)-C_gca(1,2)).^2);
                                    [~, I] = min(D); % return the # of the bead rejected
                                    
                                    [~,trajectories] = Rename_Beads(I,NROIsTotal,trajectories); % Remove the trajectories of the discarded bead
                                    [~,trajectories_sg] = Rename_Beads(I,NROIsTotal,trajectories_sg); % Same thing for the smoothed trajectory
                                    LBead(I) = [];
                                    [NROIsTotal, Beads] = Rename_Beads(I, NROIsTotal, Beads); % Remove the field corresponding to the discarded bead
                                    
                                    ZeroX_all(I) = [];
                                    ZeroY_all(I) = [];
                                    
                                    close % Close the current figure
                                    break
                                else
                                    warndlg('On the image I said!')
                                    uiwait
                                end
                                
                            case 'alt'
                                Stop_Removing_Bead = questdlg('Are you done removing beads?','Done?','Yes','No','Yes');
                                switch Stop_Removing_Bead
                                    case 'Yes'
                                        Proceed = 1;
                                        break
                                end
                        end
                    end
                end
                
            case 'No'
                Proceed = 1;
        end
        
    end
    
    %% Recalculate the origin of the curves used for the drift (sometimes, beads that are discarted were creating
    %% a shift). This part is only used to make sure that at T=0s, the curves are starting at zero in average, therefore
    %% not introducing an offset in the position of the fluorescent events detected.
    
    ZeroX = mean(ZeroX_all);
    ZeroY = mean(ZeroY_all);
    for n = 1 : NROIsTotal
        NameField = strcat('Bead_',num2str(n));
        T_sg = trajectories_sg.(NameField);
        
        frames = T_sg{1,1};
        X = T_sg{1,2} - ZeroX;
        Y = T_sg{1,3} - ZeroY;
        
        trajectories_sg.(NameField) = {frames, X, Y};
    end
    
    %% In the last part, all the curves are average together in order to calculate a proper reference curve
    
    if NROIsTotal == 0
        warndlg('There is no bead left')
        uiwait
        return
    end
    
    RefX = zeros(1,MaxFrame);
    RefY = zeros(1,MaxFrame);
    
    for n = 1 : NROIsTotal
        NameField = strcat('Bead_',num2str(n));
        T_sg = trajectories_sg.(NameField);
        RefX = RefX + T_sg{1,2};
        RefY = RefY + T_sg{1,3};
    end
    
    RefX = RefX/NROIsTotal;
    RefY = RefY/NROIsTotal;
    
    %Here we are plotting the trajectories of all the other beads after
    %subtracting the positions of the reference bead
    sigmaX = zeros(1,NROIsTotal);
    sigmaY = zeros(1,NROIsTotal);
    
    axes(h.axes2), hold off
    
    for n = 1 : NROIsTotal
        hold on
        NameField = strcat('Bead_',num2str(n));
        T_sg = trajectories_sg.(NameField);
        
        if NROIsTotal~=1
            
            if n > size(Colors,1)
                k = n - floor(n/size(Colors,1))*size(Colors,1);
            else
                k = n;
            end
            
            X_sg = T_sg{1,2} - RefX;
            Y_sg = T_sg{1,3} - RefY;
            plot(X_sg,Y_sg,'Color',Colors(k,:))
            sigmaX(n) = std(X_sg)*q;
            sigmaY(n) = std(Y_sg)*q;
        else
            
            X_sg = T_sg{1,2};
            Y_sg = T_sg{1,3};
            plot(X_sg,Y_sg,'Color',Colors(1,:))
            sigmaX = std(X_sg)*q;
            sigmaY = std(Y_sg)*q;
        end
    end
    xlabel('corrected drift along x')
    ylabel('corrected drift along y')
    title(strcat('sigmaX=',num2str(mean(sigmaX)),'nm and sigmaY=', num2str(mean(sigmaY)),'nm'));
    axis image
    Ref_position = [RefX', RefY'];
else
    warndlg('There is no beads to correct the drift with!','!!Warning!!');
    Ref_position =[];
end

h.Ref_position = Ref_position;




