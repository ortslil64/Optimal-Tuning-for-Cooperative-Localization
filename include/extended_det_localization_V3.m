function [ alpha ] = extended_det_localization_V3( PF_array,n )

for ii = 1:3
   Tc(:,:,ii) = inv(cov(PF_array.X{ii,n})); 
end
P = sum(Tc,3);
Ts = 0;
for ii = 1:3
   Ts = Ts + det(Tc(:,:,ii)) - det(P- Tc(:,:,ii));
end
for ii =1:3
   alpha(ii) = (det(P)+det(Tc(:,:,ii))-det(P-Tc(:,:,ii)))/(3*det(P) + Ts);
end



end

