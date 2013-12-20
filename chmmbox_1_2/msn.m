
num=1219
name='msn'


for i=1:num
   for j=1:num
      pinggugu(i,j)=-500;
   end
end

for i=1:num
 
    for j=1:num
        
        filename=strcat('flows/fenge/',name,'/',num2str(i),'.txt')
        fid=fopen(filename,'r');
        a=fscanf(fid,'(%f,%f)',[2,inf]);
        xunlian=a';
        fclose(fid);
        filename=strcat('flows/fenge/',name,'/',num2str(j),'.txt')
        fid=fopen(filename,'r');
        b=fscanf(fid,'(%f,%f)',[2,inf]);
        shibie=b';
        fclose(fid);   
        try
            pinggugu(i,j)=tv(xunlian,shibie);
        catch
        end
            
    end 
end


save msn.mat pinggugu

