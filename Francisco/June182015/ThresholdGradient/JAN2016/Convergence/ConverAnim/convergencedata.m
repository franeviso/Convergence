CS;
converC;
for h = 1:100;
q = linspace(0,1,100);
%q = [0.8 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.9 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99 0.995];
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
