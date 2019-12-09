function Pa = Pv2Pa(Pv,r,dr)

Pa = zeros(numel(r),1);
for i = 1:numel(r)-1
    Pa(i) = sum(r(i)./sqrt(r(i+1:end).^2-r(i)^2).*Pv(i+1:end));
end
Pa = Pa/sum(Pa)/dr;

end