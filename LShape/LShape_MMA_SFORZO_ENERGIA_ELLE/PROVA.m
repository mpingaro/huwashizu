A = [1 3 5]
B = 1:20;

somma = 0;

for i=B
    if ismember(i,A)
        i
    end
end

