close all;

origNums = 1:length(cost);
[cost, I] = sort(cost);
origNums = flip(origNums(I));
cost = flip(cost);
refs = flip(refs(I));
scatter(1:numAlgs,cost)
%scatter(1:numAlgs,cost,'.')
set(gca,'YScale','log')

cellstr(values(transMap,num2cell(transpose(intersect(cell2mat(refs(1)),cell2mat(refs(2)))))))
cellstr(values(transMap,num2cell(transpose(intersect(cell2mat(refs(numAlgs-1)),cell2mat(refs(numAlgs)))))))

%print -r100 -depsc tmp.eps;
