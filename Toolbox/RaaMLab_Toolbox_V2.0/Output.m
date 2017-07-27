function []=Output(filename,matrix)   
[n,m]=size(matrix);
fp = fopen(['../Results/',filename,'.txt'],'wt');
for i=1:n
 fprintf(fp, '%f\t', matrix(i,:));
 fprintf(fp,'\n');
end
fclose(fp);