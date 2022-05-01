function [ ary2 ] = Delete( ary1, idx )
%   delete ary1(idx)
%   return ary2 without ary1(idx)
  
    ary2 = zeros(1,length(ary1)-1);
     
    ary2(1:idx-1) = ary1(1:idx-1);
    ary2(idx:end) = ary1(idx+1:end);
  
end