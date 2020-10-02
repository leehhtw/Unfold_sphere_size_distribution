function J = costfunction2(Pa0,Pv,r,dr,weights)

Pv = Pv/sum(Pv)/dr;
Pa = Pv2Pa(Pv,r,dr);
J = (Pa-Pa0).*weights;

end