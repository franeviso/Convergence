CS = zeros(24,24);
C = zeros(24,24);
TraitAnimalData1;%trait 24x2 1 ID , 2 trait
PhylogenySorteddna100bp;%DNA-based matrix
M = 1 - M;%calculate similarity
M(logical(eye(size(M)))) = 0;
%PhylogenySorted;%DNA-based similarity matrix
%cophenetic similarity
%for i = 1:24;
%    for j = 1:24;
%        if i == j;
%           M(i,j) == 1;
%        else
%           M(i,j) = (100 - M(i,j))/100;
%        end
%    end
%end
%DNA-similarity matrix: leave as it is

%convergence matrix
for i = 1:23;
    for j = i+1:24;
        if i == j;
           C(i,j) == 1;
        else
           C(i,j) = (100 - abs(T(i,1) - T(j,1)))/100;%convergence 1
           %C(i,j) =  ((min(T(i,1),T(j,1))*100)/max(T(i,1),T(j,1)))/100;%convergence 2
           C(j,i) = C(i,j);
        end
    end
end
countNS = 0;converC=0;
%convergence events
for c = 1:24;
    sister = find(M(c,1:24) == max(M(c,1:24)));
    %M(c,1:24)%M(c,sister)%c%pause
    for d = 1:24;
%c
%d
%sister
%M(c,d)
%M(c,sister)
%pause
        if M(c,d) == M(c,sister);
           CS(c,d) = 0;CS(d,c) = 0;
        else
           if C(c,d) > C(c,sister);%convergence event 
%C(c,d)
%C(c,sister)
%pause
              CS(c,d) = C(c,d);
              CS(d,c) = C(c,d);
           else
             CS(c,d) = 0;CS(d,c) = 0;%sister closer
             countNS = countNS + 1;
             converC = countNS;
           end 
        end
   end
end

