


clear all
hold=0;


num=686;
name='aom';
tpye=1;


for j=1:num
    filename=strcat('flows/fenge/aom/170.txt')
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
        trainoutput(1,j)=tv(xunlian,shibie);
    catch
         trainoutput(1,j)=hold;
    end
    
    
    filename=strcat('flows/fenge/edonkey/879.txt')
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
        trainoutput(2,j)=tv(xunlian,shibie);
    catch
         trainoutput(2,j)=hold;
    end   
 
    
    filename=strcat('flows/fenge/http/20.txt')
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
        trainoutput(3,j)=tv(xunlian,shibie);
    catch
         trainoutput(3,j)=hold;
    end   
    
    filename=strcat('flows/fenge/msn/515.txt')
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
        trainoutput(4,j)=tv(xunlian,shibie);
    catch  
         trainoutput(4,j)=hold;
    end   
       
    filename=strcat('flows/fenge/pplive/118.txt')
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
        trainoutput(5,j)=tv(xunlian,shibie);
    catch 
         trainoutput(5,j)=hold;
    end    
    
    filename=strcat('flows/fenge/smtp/630.txt')
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
        trainoutput(6,j)=tv(xunlian,shibie);
    catch 
         trainoutput(6,j)=hold;
        
    end    
    
    
    filename=strcat('flows/fenge/vedio/487.txt')
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
        trainoutput(7,j)=tv(xunlian,shibie);
    catch 
        trainoutput(7,j)=hold;
    end
    
    
end 

[r,c]=max(trainoutput);

right=0;
for i=1:num
    if c(i)==type
        right=right+1;
    end
end
zhunql(type)=right/num


%clear r
%clear c


