for r = 1:10;
plantcompLast;
animcompLast;

m=69;indX = randperm(round(nnz(CP(r,:))),m);
n=24;indY = randperm(round(nnz(CA(r,:))),n);
CP = indX;CA = indY;

q = linspace(0,1,100);
for h = 1:100;
count = 0;C = zeros(1,1);%Number of replicates!!
    for i = 1:length(CA);
        for j = 1:length(CP);
            if CA(1,i) > 0 & CP(1,j) > 0;
               count = count + 1;
%CA(r,i)
%CP(r,j)
%q(1,h)
               C(count,1) = ((min(CA(1,i),CP(1,j))*100)/max(CA(1,i),CP(1,j)))/100;
%pause
               %C(count,1) = Wabs(CA(r,i) - CP(r,j));
               %C(count,1) = (100 - abs(CA(r,i) - CP(r,j)))/100;
            end
        end
    end
   Events = find(C(:,1) > q(1,h));
   PrC(h,1) = length(Events)/count;
%pause
%end
%comp = mean(PrC);
%mincomp = min(PrC);
%maxcomp = max(PrC);


%fid = fopen('meancomplementariedad.txt','a');fprintf(fid,'%6.8f %6.8f\n',q(1,h),comp);fclose(fid); 
%fid = fopen('mincomplementariedad.txt','a');fprintf(fid,'%6.8f %6.8f\n',q(1,h),mincomp);fclose(fid); 
%fid = fopen('maxcomplementariedad.txt','a');fprintf(fid,'%6.8f %6.8f\n',q(1,h),maxcomp);fclose(fid); 
end
fid = fopen('complementarityAll.txt','a');fprintf(fid, [repmat('% 6f',1,size(PrC(:,1),1)), '\n'],PrC);fclose(fid);
end%r
