% function for stochastic matrix dP
function dy = ss_dy(t,y,Pdim1,Pdim2)

global r_ std Qmat odiag;

P = reshape(y,Pdim1,Pdim2);

% allow to mutate to NNs with probs odiag
for i=2:Pdim1-1 % leave out edges...
    for j=2:Pdim2-1
            odiag = mod(abs(std*randn(4,1)),0.25);
            Qmat(i,j,i+1,j) = odiag(1);
            Qmat(i,j,i-1,j) = odiag(2);
            Qmat(i,j,i,j+1) = odiag(3);
            Qmat(i,j,i,j-1) = odiag(4);
        Qmat(i,j,i,j) = 1 - sum(odiag);
    end
end

% calculate dP according to mut matrix eqn
dP = zeros(Pdim1,Pdim2);
for i=1:Pdim1
    for j=1:Pdim2
        dP(i,j) = r_.*squeeze(sum(sum(P.*squeeze(Qmat(:,:,i,j)))));
    end
end

total = sum(sum(Qmat(4,3,:,:),3),4);
total

t

dy = reshape(dP,Pdim1*Pdim2,1);
