function [J, dJ] = costfunction(Pa0,Pv,r,dr,lambda,weights)

Pv = Pv/sum(Pv)/dr;
Pa = Pv2Pa(Pv,r,dr);
J = 1/2/numel(r)*sum((Pa-Pa0).^2.*weights.^2) + 1/2/numel(r)*lambda*sum(Pv.^2);

dJ = zeros(numel(Pv),1);
eps = 1e-10;
for i = 1:numel(Pv)
    Pv1 = Pv; Pv2 = Pv;
    Pv1(i) = Pv(i) + eps;
    Pv2(i) = Pv(i) - eps;
    Pa1 = Pv2Pa(Pv1,r,dr);
    Pa2 = Pv2Pa(Pv2,r,dr);
    dPa = (Pa1-Pa2)/2/eps;
    dJ(i) = 1/numel(r)*sum((Pa-Pa0).*dPa.*weights.^2);
end
dJ = dJ + 1/numel(r)*lambda*Pv;

end