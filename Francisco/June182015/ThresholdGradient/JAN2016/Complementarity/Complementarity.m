%First test outputs
%AnimComp;
%PlantComp;
%Second test outputs
plantcompLast;
animcompLast;

%# pick random columns plants 69
m=69;
indX = randperm( size(CP,1) );
indX = indX(1:m);

%# pick random columns anim 24
n=24;
indY = randperm( size(CA,1) );
indY = indY(1:n);

%# filter data
CP = CP(1,indX)

CA = CA(1,indY)

size(CP)
size(CA)
pause

q = linspace(0,1,100);
for h = 1:100;
    %q = [0.8 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.9 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99 0.995];
for r = 1:10;count = 0;C = zeros(1,1);%Number of replicates!!
    for i = 1:length(CA);
        for j = 1:length(CP);
            if CA(r,i) > 0 & CP(r,j) > 0;
               count = count + 1;
%CA(r,i)
%CP(r,j)
%q(1,h)
               C(count,1) = ((min(CA(r,i),CP(r,j))*100)/max(CA(r,i),CP(r,j)))/100;
%pause
               %C(count,1) = Wabs(CA(r,i) - CP(r,j));
               %C(count,1) = (100 - abs(CA(r,i) - CP(r,j)))/100;
            end
        end
    end
 Events = find(C(:,1) > q(1,h));
    PrC(r,1) = length(Events)/count;
%pause
end
comp = mean(PrC);
mincomp = min(PrC);
maxcomp = max(PrC);
fid = fopen('meancomplementariedad.txt','a');fprintf(fid,'%6.8f %6.8f\n',q(1,h),comp);fclose(fid); 
fid = fopen('mincomplementariedad.txt','a');fprintf(fid,'%6.8f %6.8f\n',q(1,h),mincomp);fclose(fid); 
fid = fopen('maxcomplementariedad.txt','a');fprintf(fid,'%6.8f %6.8f\n',q(1,h),maxcomp);fclose(fid); 
end
