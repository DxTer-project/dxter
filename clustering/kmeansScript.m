k = 30;

ids = kmeans(cost,k);
means = zeros(k,1);

%close all;
figure(2)

intersections = cell([k 1]);

for i=1:k
	ind = find(ids == i);
means(i) = mean(cost(ind));
        intersection = cell2mat(refs(ind(1)));
for j=2:length(ind)
	intersection=intersect(intersection,cell2mat(refs(ind(j))));
end
intersections(i) = mat2cell(transpose(intersection),[length(intersection)],[1]);
%means(i)
%cellstr(values(transMap,num2cell(transpose(intersection))))
end

[means,I] = sort(means);
means = flip(means);
intersections = flip(intersections(I));

scatter(1:k, means); hold on;
set(gca,'YScale','log')
for i=1:k
text(i, means(i), num2str(i)); hold on;
disp(['Group ' num2str(i)]);
disp(['         Contains ' num2str(length(find(ids==i)))]);
if (i > 1)
  disp(['Improvement of the following group over previous = ' num2str(1-means(i)/means(i-1))]);
end
cellstr(values(transMap,num2cell(cell2mat(intersections(i)))))
end

maxGroupID = 1:k;
maxGroupID = flip(maxGroupID(I));
maxGroupID = maxGroupID(k);
lastGroup = find(ids == maxGroupID);
id=fopen('transformations.txt','w');
fprintf(id,'printing %d groups of transformations for last group\n',length(lastGroup));
for i=1:length(lastGroup)
	fprintf(id,'%d:\n',origNums(lastGroup(i)));
fprintf(id,'    cost %e\n',cost(lastGroup(i)));
	tmp = cellstr(values(transMap,num2cell(cell2mat(refs(lastGroup(i))))));
	for j=1:length(tmp)
		tmp2 = tmp(j);
		tmp2 = tmp2{1};
                fprintf(id, '%s\n',tmp2);
        end
	fprintf(id, '\n\n');
end
fclose(id);

c= [1 1 0;
      1 0 1;
      0 1 1;
      1 0 0;
      0 1 0;
      0 0 1;
      0 0 0
    .5 .5 0
    0 .5 .5
    .5 0 .5
    .5 0 0
    0 .5 0
    0 0 .5
    .5 .5 .5];

while (k > length(c(:,1))) 
  disp 'Reusing colors';
  c = [c;c];
end

figure(3);

scatter(1:numAlgs,cost,4,c(ids,:));
set(gca,'YScale','log')


