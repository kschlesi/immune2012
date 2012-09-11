
function dy = Pmutn(t,y,Pdim1,Pdim2)

global r_ mrate ijs;

P = reshape(y,Pdim1,Pdim2);

% random walk steps of each location, added to location itself, for matrix
% of new locations
rwsteps = cast(unifrnd(-1.5,1.5,Pdim1,Pdim2,2),'int16') + ijs;

Pmults = zeros(Pdim1,Pdim2); % scalar field of #pathogen mutating INTO each P site
for i=0:Pdim1-1
    for j=0:Pdim2-1
        if squeeze(rwsteps(mod(i+1,Pdim1)+1,j+1,:))==[i+1;j+1]
            Pmults(i+1,j+1) = Pmults(i+1,j+1)+P(mod(i+1,Pdim1)+1,j+1);
        end
        if squeeze(rwsteps(mod(i+Pdim1-1,Pdim1)+1,j+1,:))==[i+1;j+1]
            Pmults(i+1,j+1) = Pmults(i+1,j+1)+P(mod(i+Pdim1-1,Pdim1)+1,j+1);
        end
        if squeeze(rwsteps(i+1,mod(j+1,Pdim2)+1,:))==[i+1;j+1]
            Pmults(i+1,j+1) = Pmults(i+1,j+1)+P(i+1,mod(j+1,Pdim2)+1);
        end
        if squeeze(rwsteps(i+1,mod(j+Pdim2-1,Pdim2)+1,:))==[i+1;j+1]
            Pmults(i+1,j+1) = Pmults(i+1,j+1)+P(i+1,mod(j+Pdim2-1,Pdim2)+1);
        end
    end
end

t

dP = mrate.*(Pmults-P)+r_.*P;
dy = reshape(dP,Pdim1*Pdim2,1);
