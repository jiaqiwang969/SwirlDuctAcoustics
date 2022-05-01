function [G_nm]=cheb_cumKxCell_r(Gr_nm,Ratio)
%initialization
G_nm=chebfun(0,[Ratio, 1]);
for kkk=1:length(Gr_nm) 
    G_nm=G_nm+Gr_nm{kkk};
end
