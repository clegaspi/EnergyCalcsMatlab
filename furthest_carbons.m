function idx = furthest_carbons(path_to_file)
[xyz, atom_type]=ampac_to_xyz(path_to_file);
c_idx=find(cellfun(@(x)strcmpi(x,'C'),atom_type));
nc=length(c_idx);
dist=zeros(nc);
vectsum = zeros(1,3);
centroid=vectsum;
for i = 1:nc
    for j = i+1:nc
        dist(i,j)= sum((xyz(c_idx(i),:)-xyz(c_idx(j),:)).^2)^0.5;
        vectsum = vectsum + (xyz(c_idx(i),:)-xyz(c_idx(j),:))./norm(xyz(c_idx(i),:)-xyz(c_idx(j),:));
    end
    centroid = centroid + xyz(c_idx(i),:);
end
[~,i]=max(dist(:));
[i,j]=ind2sub(size(dist),i);

centroid = centroid ./ nc;
carb_vect = xyz(c_idx(i),:)-xyz(c_idx(j),:);
carb_vect = carb_vect ./ norm(carb_vect);
vectsum = vectsum ./ norm(vectsum);
err_diff = (vectsum-carb_vect)./vectsum;

disp(['The furthest carbons are ',num2str(c_idx(i)),' and ',num2str(c_idx(j))]);
disp(['The vector formed by these carbons is ( ',...
    num2str(carb_vect(1)),' , ',...
    num2str(carb_vect(2)),' , ',...
    num2str(carb_vect(3)),' )']);
disp(['The normalized vector sum is ( ',...
    num2str(vectsum(1)),' , ',...
    num2str(vectsum(2)),' , ',...
    num2str(vectsum(3)),' )']);
disp(['The error difference between the two vectors is ', num2str(err_diff*100),'%']);
idx = [c_idx(i), c_idx(j)];
disp(['The centroid is at ( ',num2str(centroid),' )']);
end