function [A] = binread(fname)
%BINREAD Summary of this function goes here
%   Detailed explanation goes here
fid=fopen(fname,'rb');
N=fread(fid,1,'int','l');
M=fread(fid,1,'int','l');
A=zeros(N,M);
for j=1:M
   A(:,j)=fread(fid,N,'double','l');
end
fclose(fid);




end

