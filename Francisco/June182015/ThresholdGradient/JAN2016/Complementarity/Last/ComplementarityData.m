TraitAnimalData1;
TraitPlantData1;
for h = 1:100;
q = linspace(0,1,100);
    %q = [0.75 0.77 0.8 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.9 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99 0.995];
%for r = 1:15;
count = 0;C = zeros(1,1);Events = 0;
    for i = 1:length(CA);
        for j = 1:length(CP);
            if CA(i,1) > 0 & CP(j,1) > 0;
               count = count + 1;
               C(count,1) = ((min(CA(i,1),CP(j,1))*100)/max(CA(i,1),CP(j,1)))/100;
               %C(count,1) = (100 - abs(CA(i,1) - CP(j,1)))/100;
            end
        end
    end
    %C(:,1)/max(C)
    Events = find(C(:,1)/max(C) > q(1,h));
%q(1,h)
    PrC = length(Events)/count;
%end
%comp = mean(PrC);
%mincomp = min(PrC);
%maxcomp = max(PrC);
fid = fopen('meancomplementariedadData1.txt','a');fprintf(fid,'%6.8f %6.8f\n',q(1,h),PrC);fclose(fid); 
%fid = fopen('mincomplementariedad.txt','a');fprintf(fid,'%6.8f %6.8f\n',q(1,h),mincomp);fclose(fid); 
%fid = fopen('maxcomplementariedad.txt','a');fprintf(fid,'%6.8f %6.8f\n',q(1,h),maxcomp);fclose(fid); 
end