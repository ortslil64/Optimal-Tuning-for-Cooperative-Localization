function [ alpha ] = det_localization_V3( PF_array,n )

for ii = 1:3
   Tc(ii) = 1/det(cov(PF_array.X{ii,n})); 
end
for ii =1:3
   alpha(ii) = Tc(ii)./sum(Tc);
end



end

