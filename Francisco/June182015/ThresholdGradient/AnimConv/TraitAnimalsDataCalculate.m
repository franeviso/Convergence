CS = zeros(40,40);
fileID = fopen('TraitAnimalData2.csv');
in = fscanf(fileID,'%c');fclose(fileID);
A = regexp(in,'(^,|)([^\n]+)', 'match');
w1 = char(A);count = 0;count1=0;
for j = 1:length(A);
    count = count + 1;
    wHP = regexp(w1(j,:),' ');
    wordP = w1(j,1:wHP(1,1));
    cP(count,:) = cellstr(wordP);
end

for j = 1:length(A);
    count1 = count1 + 1;
    wHC = regexp(w1(j,:),',');
    numberP = w1(j,wHC(1,1)+1:wHC(1,1)+4);
    cN(count1,:) = cellstr(numberP);
end
matdata=cellfun(@str2num,cN);
%b = unique_no_sort(cP);
for i = 1:length(cP)-1;ind = zeros(1,1);
    ind = find(ismember(cP,cP(i)));
    for j = i+1:length(cP);
        if isempty(ind);
matdata(i,1)
matdata(j,1)
kk = ((min(matdata(i,1),matdata(j,1))*100)/max(matdata(i,1),matdata(j,1)))/100
pause

           %CS(i,j) = (100 - abs(matdata(i,1) - matdata(j,1)))/100;
           CS(i,j) = ((min(matdata(i,1),matdata(j,1))*100)/max(matdata(i,1),matdata(j,1)))/100;
           CS(j,i) = CS(i,j);
        else
           %trait value nosister closer to focal than sister
           ind(ind==i) = [];
               for k = 1:length(ind);
                   %if ind(k,1) ~= j;
                      if abs(matdata(ind(k,1),1) - matdata(i,1)) < abs(matdata(j,1) - matdata(i,1));
%i
%j
%ind(k,1)
%matdata(ind(k,1),1)
%matdata(i,1)
%matdata(j,1)
%pause
                   CS(i,j) = 0;
                   CS(j,i) = 0;
                   else
                   CS(i,j) = ((min(matdata(i,1),matdata(j,1))*100)/max(matdata(i,1),matdata(j,1)))/100;
                   %CS(i,j) = (100 - abs(matdata(i,1) - matdata(j,1)))/100;
                   CS(j,i) = CS(i,j);
                   end
                   %end
               end   
        %sister = find(ind == j)
        %if isempty(sister);
        %   CS(i,j) = (100 - abs(matdata(i,1) - matdata(j,1)))/100;
        %   CS(j,i) = CS(i,j);
        %else
          %if distance nosister-focal larger than sister-focal
        %   sisterfocal = find(
        end
    end 
end

