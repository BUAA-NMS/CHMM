function [ pinggu ] = tv( xunlian , shibie )

Z=xunlian;
Z1=shibie;

%clear all
load chmmsim
%T=length(data.Xseries);
%T=132;









s=size(Z);
T=s(1:1,end-1);
X=Z(1:T,end-1);
Y=Z(1:T,end);




s1=size(Z1);
T1=s1(1:1,end-1);
X1=Z1(1:T1,end-1);
Y1=Z1(1:T1,end);




% Train up GMM on this data
chmm=struct('hmmone',[],'hmmtwo',[]);
chmm.hmmone.K=2; 
chmm.hmmtwo.K=2;

%chmm=chmmrandinit(X,Y,chmm,['full';'full']);
chmm=chmminit(X,Y,chmm,['full';'full']);
% $$$ chmm.hmmone.state(1).Mu=simchmm.hmmone.state(2).Mu';

%disp('Means of CHMM initialisation');
chmm.hmmone.state(:).Mu;
chmm.hmmtwo.state(:).Mu;

% Train up HMM on observation sequence data using Baum-Welch
% This uses the forward-backward method as a sub-routine
%disp('We will now train the CHMM using Baum/Welch');
%disp(' ');
%disp('Press a key to continue');
%pause
%disp('Estimating CHMM');

chmm.train.cyc=20;
chmm.hmmone.obsmodel='Gauss';
chmm.hmmtwo.obsmodel='Gauss';
chmm.hmmone.train.obsupdate=ones(1,chmm.hmmone.K);% update observation models ?
chmm.hmmtwo.train.obsupdate=ones(1,chmm.hmmtwo.K);% update observation models ?
chmm.hmmone.train.init=1;         % Yes, we've already done initialisation
chmm.hmmtwo.train.init=1;         % Yes, we've already done initialisation


chmm=chmmtrain(X,Y,T,chmm);
%disp('Means');
%chmm.hmmone.state.Mu
%chmm.hmmtwo.state.Mu


%disp('Initial State Probabilities, Pi');
chmm.Pi;
%disp('State Transition Matrix, P');
chmm.P;



%disp('We will estimate the Viterbi path of the CHMM, ');
%disp(' ');
%disp('Press a key to continue');
%pause

chmm.P=chmm.P';                        % Will's code uses transpose 
% joint viterbi path
%%%%%%%%%%%%%%%%%%%%%%%%%以上是模型训练部分%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp('模型训练成功');


[block.joint]=chmmdecode(X1,Y1,T1,chmm);
%[block.joint]=chmmdecode(X,Y,T,chmm);

% We are often interested in the Viterbi path of the individual chains. We can 
% estimate these in two ways: a) marginalise the transition probabilities and 
% run decoding for each of the two HMMs; b) marginalise the joint state
% sequence directly 


block.jointviterbi=block.joint.q_star;
block.hmmoneviterbi=zeros(size(block.jointviterbi));
block.hmmtwoviterbi=zeros(size(block.jointviterbi));
clvlabels=[chmm.hmmone.K 1:chmm.hmmone.K-1];
clv=mod(block.joint.q_star,chmm.hmmone.K);
for i=0:max(clv),
  ndx=find(clv==i);
  block.hmmoneviterbi(ndx)=clvlabels(i+1)+zeros(size(ndx));
end;
clvlabels=[chmm.hmmtwo.K 1:chmm.hmmtwo.K-1];
clv=mod(ceil(block.joint.q_star/chmm.hmmone.K),chmm.hmmtwo.K);
for i=0:max(clv),
  ndx=find(clv==i);
  block.hmmtwoviterbi(ndx)=clvlabels(i+1)+zeros(size(ndx));
end;
clear ndx clv;
  



% single chain viterbi path
  % priors for individual chains
  Pone=ones(1:chmm.hmmone.K);
  Pone=Pone./length(Pone);
  Ptwo=ones(1:chmm.hmmtwo.K);
  Ptwo=Ptwo./length(Ptwo);

  % first chain; integrate out 2nd chain
  hmm=chmm.hmmone;
  hmm.P=mdsum(mdprod(hmm.P,Ptwo,3),3)';
  [block.hmmone]=hmmdecode(X1,T1,hmm);

  % second chain; integrate out 1st chain
  hmm=chmm.hmmtwo;
  hmm.P=mdsum(mdprod(hmm.P,Pone,2),2)';
  [block.hmmtwo]=hmmdecode(Y1,T1,hmm);
   
  
  block.joint.q_star;
  
  
  
  % 评估值是log和
  pinggu=0;
  for i=[1:T1-1]
  pinggu=pinggu+log(chmm.P(block.joint.q_star(i),block.joint.q_star(i+1)));
  end
  pinggu 

  
  
  %___________________________________________
  
  
    
% Z2=[53.982,10;53.982,10;53.981,40;21.335,10;53.979,10;53.982,10;53.981,10;53.981,40;21.335,10;53.979,10;53.981,10;53.981,10;53.981,40;22.405,10;53.980,24;21.173,10;53.979,24;21.239,10;53.978,24;23.424,10;46.977,40;53.012,24;21.239,10;46.982,24;20.212,24;18.633,24;42.186,15;20.969,57;22.553,51;23.560,16;19.395,19;19.494,37;34.873,19;21.239,40;24.728,19;22.095,134;33.941,19;20.792,109;27.679,19;20.531,43;19.294,19;19.085,43;19.085,19;19.085,43;19.031,19;19.085,43;19.731,19;20.899,43;19.191,19;19.031,43;19.085,19;19.085,43;19.085,19;19.031,43;19.031,19;20.043,43;19.138,19;19.685,43;20.645,19;19.912,43;19.085,19;19.294,43;19.138,19;21.106,43;19.912,19;19.294,43;19.191,19;19.542,43;19.685,19;20.492,43;19.294,19;19.345,25;19.243,19;19.345,46;21.038,27;29.703,12;43.447,12;51.744,10;53.981,10;53.982,10;24.928,16;53.975,10;49.217,12;42.241,12;51.757,10;53.982,10;42.251,16;53.679,10;53.981,10;50.010,12;42.229,12;51.245,10;45.244,16;53.357,10;53.981,10;53.981,10;45.250,16;49.207,12;24.786,12;51.239,10;53.982,10;53.981,10;47.010,16;53.008,10;50.674,12;42.241,12;50.665,10;53.981,10;47.012,16;53.008,10;53.981,10;51.254,12;24.757,12;50.657,10;47.004,16;53.357,10;53.982,10;53.981,10;47.005,16;48.236,12;25.453,12;51.237,10;53.682,10;53.981,10;49.217,16;52.217,10;51.765,12;24.829,12;49.992,10;53.981,10;50.006,16;52.218,10;53.682,10;51.761,12;26.551,12;50.659,10;49.213,16;51.755,10;53.982,10;53.981,10;50.677,16;46.986,12;25.416,12;49.187,10;53.982,10;54.261,10;50.007,16;51.246,10;52.225,12;42.224,12;49.201,10;53.982,10;50.675,16;51.247,10;53.981,10;52.228,12;42.197,12;48.228,10;50.770,16;51.163,10;53.981,10;53.981,10;51.258,16;46.990,12;42.220,12;46.977,10;53.981,10;53.981,10;52.225,16;49.203,10;52.638,12;45.240,12;45.196,10;53.981,10;52.225,16;49.202,10;53.981,10;53.016,12;42.232,12;45.201,10;52.638,16;48.234,10;53.981,10;53.981,10;52.638,16;42.235,12;45.213,12;42.201,10;53.981,10;53.981,10;53.016,16;46.995,10;53.362,12;42.181,12;42.216,10;53.978,10;53.363,16;45.212,10;53.982,10;53.683,12;25.224,12;42.151,10;53.682,16;42.179,10;53.982,10;53.981,10;53.685,16;42.198,12;22.577,10;24.683,12;53.972,10;53.981,10;53.981,10;25.599,16;53.976,10;25.119,12;42.171,12;53.679,10;53.981,10;42.259,16;53.679,10;53.984,10;42.255,12;42.182,12;53.357,10;42.851,16;53.633,10;54.262,10;53.982,10;42.240,12;34.776,12;41.437,16;53.353,10;53.981,10;53.982,10;45.257,16;53.355,10;45.251,12;34.704,12;53.297,10;53.982,10;47.005,16;53.008,10;53.981,10;47.005,12;34.736,12;52.944,10;48.249,16;52.631,10;53.981,10;53.981,10;48.247,12;34.654,16;21.959,12;52.559,10;53.981,10;53.981,10;49.220,16;52.216,10;48.246,12;42.828,12;52.153,10;53.981,10;50.008,16;51.759,10;53.981,10;48.248,12;45.546,12;51.686,10;50.678,16;51.246,10;53.981,10;53.981,10;49.215,12;45.558,16;22.148,12;51.159,10;53.981,10;53.982,10;51.256,16;50.664,10;50.005,12;45.560,12;51.162,10;53.981,10;51.767,16;49.997,10;53.981,10;50.005,12;45.561,12;50.570,10;51.770,16;49.990,10;53.981,10;53.982,10];

 
 
 
 