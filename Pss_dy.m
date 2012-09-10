function dP = Pss_dy(t,P,Pdim1)

global r_ D;

dmut = zeros(Pdim1,1);
% boundaries - forward differences?
for i=2:Pdim1-1
    dmut(i) = D*(P(i-1)-2*P(i)+P(i+1));
end

dP = dmut;
%dP = dmut + r_.*P;
