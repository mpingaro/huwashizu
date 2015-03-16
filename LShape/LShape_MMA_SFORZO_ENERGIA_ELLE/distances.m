function [distanze] = distances(coordinates,element)

nele = length(element);
nnod = length(coordinates);
bari = zeros(nele,2);
distanze = zeros(nele);

bari(:,1) = coordinates(element(:,1),1)+coordinates(element(:,2),1)+coordinates(element(:,3),1);
bari(:,2) = coordinates(element(:,1),2)+coordinates(element(:,2),2)+coordinates(element(:,3),2);

bari = bari./3;

for j=1:nele
    distanze(:,j) = ((bari(:,1)-repmat(bari(j,1),nele,1)).^2+(bari(:,2)-repmat(bari(j,2),nele,1)).^2).^(1/2);
end

