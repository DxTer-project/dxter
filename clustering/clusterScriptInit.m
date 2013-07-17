%y = pdist(cost);
z = linkage(cost,'ward','euclidean','savememory','on');

save('clusterInitResults');
