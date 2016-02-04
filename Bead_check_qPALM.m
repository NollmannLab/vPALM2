function bead = Bead_check_qPALM(bead, Max_Frames, AnalysisMethods, Applied3D)

bead(:,bead(1,:)>Max_Frames) = [];
for k = 1 : Max_Frames
    
     Check = find(bead(1,:)==k); %Return the position in the matrix bead of the events detected by MTT in the cropped area for the frame #i
    
    %In case there is no event detected at frame #i, an event will be automatically created by averaging over the surrounding events
    if size(Check,2) == 0
        if k<5
            X = mean(bead(2,1:5));
            Y = mean(bead(3,1:5));
            
            if (AnalysisMethods == 7 || AnalysisMethods == 5)  && Applied3D
                Z = mean(bead(4,1:5));
            end
            
        elseif k<=size(bead,2) - 2
            X = mean(bead(2,k-2:k+2));
            Y = mean(bead(3,k-2:k+2));
            
            if (AnalysisMethods == 7 || AnalysisMethods == 5) && Applied3D
                Z = mean(bead(4,k-2:k+2));
            end
            
        else
            X = mean(bead(2,k-5:k-1));
            Y = mean(bead(3,k-5:k-1));
            
            if (AnalysisMethods == 7 || AnalysisMethods == 5) && Applied3D
                Z = mean(bead(4,k-5:k-1));
            end
            
        end
        
        if (AnalysisMethods == 7 || AnalysisMethods == 5) && Applied3D
            Elementk = [k; X; Y; Z];
        else
            Elementk = [k; X; Y];
        end
        bead = [bead(:,1:k-1) , Elementk , bead(:,k:end)];
    end
end



