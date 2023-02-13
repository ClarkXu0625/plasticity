%step 1
Ncell = 60;
Group = randi(4, 1, Ncell);

%step 2
C = rand(Ncell) < 0.2;
%C = zeros(Ncell);
for i = 1:Ncell
    for j = 1:Ncell
        if Group(i) == Group(j)
            x = rand;
            if x<0.5
                C(i, j) = 1;
            end
        end
    end
end

%step 4
corrColumn = C(:, 1);
for i = 2:Ncell
    temp = C(:, i);
    corrColumn = [corrColumn temp];
end
R1 = corrcoef(corrColumn);

corrRow = C(1, :);
for i = 2:Ncell
    temp = C(i, :);
    corrRow(end+1, :) = temp;
end

R2 = corrcoef(corrRow);
total_corr = R1+R2;


histogram(total_corr(:));

%step 5
sorted = sort(total_corr(:));
threshold = sorted(length(sorted)*0.9);

%step 6

groupId = zeros(1, Ncell);

GC = 0; %group count


for i = 1:Ncell
    if groupId(i) == 0
        GC = GC+1;
        remain = zeros(1, 1);
        remain(1) = i;
        used = [];
        
        %9a
        while isempty(remain) == 0 %while remain list not empty
            temp = remain(1); %temp is the first element of remaining list
            
            if length(remain) == 1
                remain = [];
            else
                remain = remain(2:end); %delete the first element
            end

            %9b
            temp1 = total_corr(temp, :);

            for j = 1:Ncell
                if temp1(j)>=threshold
                    groupId(j) = GC; %9c
                    remain(end+1) = j; %9c
                    
                    
                end
            end

            used(end+1) = i; %9d


            %step 9e remove repeated and used elements
            if isempty(remain) == 0 %when remain is not empty
                remain = sort(remain);
                %delete all element that are in the used list
                k = 1;
                while k <= length(remain)
                    if ismember(remain(k), used)
                        remain(k) = [];
                    else
                        k = k+1;
                    end
                end
                %delete repeating elements
                k = 2;
                while k <= length(remain)
                    if remain(k) == remain(k-1)
                        remain(k) = [];
                    else
                        k = k+1;
                    end
                end
            end

        end
        
        
    end        
end

disp(groupId);
