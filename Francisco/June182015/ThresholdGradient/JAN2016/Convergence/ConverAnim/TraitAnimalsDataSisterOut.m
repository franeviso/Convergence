CS = zeros(24,24);
fileID = fopen('TraitAnimalData1_24.csv');
in = fscanf(fileID,'%c');fclose(fileID);
A = regexp(in,'(^,|)([^\n]+)', 'match');
w1 = char(A);count = 0;count1=0;countNS = 0;converC=0;
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
for i = 1:length(cP)-1;ind = zeros(1,1);
    ind = find(ismember(cP,cP(i)));
    for j = i+1:length(cP);
        if isempty(ind);%no sp same genera (compute pheno similarity)
           CS(i,j) = ((min(matdata(i,1),matdata(j,1))*100)/max(matdata(i,1),matdata(j,1)))/100;
           CS(j,i) = CS(i,j);
        else
        ind(ind==i) = [];%remove i
        sister = find(ind == j);
              if isempty(sister);%nonsister
                 for k = 1:length(ind);
                     if abs(matdata(ind(k,1),1) - matdata(i,1)) < abs(matdata(j,1) - matdata(i,1));%sister closer than nonsister-sister then CS = 0
                        CS(i,j) = 0;
                        CS(j,i) = 0;
                        countNS = countNS + 1;
                        converC = countNS;
                        break
                     end
                 end
                 CS(i,j) = ((min(matdata(i,1),matdata(j,1))*100)/max(matdata(i,1),matdata(j,1)))/100;
                 CS(j,i) = CS(i,j);
              else%sister
                 CS(i,j) = 0;
                 CS(j,i) = 0;
              end
        end
   end%j
%i
%j
%CS
%pause
end%i
