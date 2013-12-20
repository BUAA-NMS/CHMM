function [ output_args ] = autotv( Z,Z1 )
%UNTITLED1 Summary of this function goes here
%  Detailed explanation goes here



fid=fopen('flows/Z.txt','r');
a=fscanf(fid,'(%f,%f)',[2,inf]);
Z=a';
fclose(fid);


fid=fopen('flows/Z1.txt','r');
b=fscanf(fid,'(%f,%f)',[2,inf]);
Z1=b';
fclose(fid);




