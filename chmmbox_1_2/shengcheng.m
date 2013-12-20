
comp0=0.1174; miu0=7.994;sigma0=0.236;
comp1=0.2135;miu1=21.71;sigma1=47.04;
comp2=0.1255;miu2=38.57;sigma2=6.561;
comp3=0.2867;miu3=47.96;sigma3=5.359;
comp4=0.2570;miu4=53.18;sigma4=1.087;


for i=1:2
    cur=rand;
    if cur>=0 && cur< comp0
        normrnd(miu0,sigma0,1,1)   
    end
    if cur>=comp0 && cur< comp0+comp1
        normrnd(miu1,sigma1,1,1)  
    end
    if cur>=comp0+comp1 && cur< comp0+comp1+comp2
        normrnd(miu2,sigma2,1,1)    
    end
    if cur>=comp0+comp1+comp2 && cur< comp0+comp1+comp2+comp3
        normrnd(miu3,sigma3,1,1)    
    end
    if cur>=comp0+comp1+comp2+comp3 && cur<=1
        normrnd(miu4,sigma4,1,1)    
    end
end

    

