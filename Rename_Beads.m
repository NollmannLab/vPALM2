function [NBeads, Beads] = Rename_Beads(Number, NBeads, Beads)

%% This function is used during the drift correction to remove from 
%% the structures "Beads", "trajectories" and "trajectoriesSG" the
%% field related to "Number". It returns a new structure with one field less
%% but the all the name going from 1 to NBeads in ascending order.

for k = Number : NBeads
	ThisField = strcat('Bead_',num2str(k));
	if k < NBeads
        NextField = strcat('Bead_',num2str(k+1));
        NextBeadValue = getfield(Beads, NextField);
        Beads = setfield(Beads, ThisField, NextBeadValue);
	elseif k == NBeads
        Beads = rmfield(Beads,ThisField);
        NBeads = NBeads -1 ;
	end
end
