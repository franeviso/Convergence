CS;
converC;
for h = 1:100;
q = linspace(0,1,100);
count2 = 0;events = 0;converA = 0;converB=0;
    for i = 1:length(CS)-1;
        for j = 1+i:length(CS);
            if CS(i,j) > 0;
               if CS(i,j) > q(1,h);
                  events = events + 1;
                  converA = events;
               else
                  count2 = count2 +1;
                  converB=count2;
               end
            end
        end
    end
    PrC(1,h) = converA/(converA+converB+converC);
fid = fopen('convergencefromdata1.txt','a');fprintf(fid,'%6.8f %6.8f\n',q(1,h),PrC(1,h));fclose(fid); 
end
