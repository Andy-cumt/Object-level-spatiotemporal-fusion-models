%%%Robert

function robertsNum = Robert(image,i)

[m,n]=size(image);
edge = zeros(m,n);
robertsNum = zeros(m,n);

for j=1:m-1
    for k=1:n-1
        robertsNum(j,k) = abs(image(j,k)-image(j+1,k+1)) + abs(image(j+1,k)-image(j,k+1));
    end
end

B = sort(robertsNum(:),'descend');
T = B(round(i * m * n));
for j=1:m
    for k=1:n
        if robertsNum(j,k) >= T
            edge(j,k) = 1;
        end
    end
end


