function [] = binwrite(fname,A)
%BINREAD Output Matrix as Raw Binary
%
fid=fopen(fname,'wb');
[M,N]=size(A);
fwrite(fid,M,'int');
fwrite(fid,N,'int');
for j=1:N
   fwrite(fid,A(:,j),'double');
end
fclose(fid);




end

