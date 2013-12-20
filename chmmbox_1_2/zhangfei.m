% A demonstration of the CHMM software using a Gaussian observation
% model on AR features
clear all
load chmmsim
%T=length(data.Xseries);
T=132;
%X=data.Xseries(1:T,:);Y=data.Yseries(1:T,:);
Z=[21.563025,1380;21.051500,1384;20.913617,1378;54.087780,1384;22.334685,201;21.051500,1378;21.440680,1380;53.960815,1384;37.832698,1384;21.797836,1378;21.440680,1380;54.055552,201;22.627578,1380;21.314789,1384;21.563025,1382;54.073929,1384;22.434527,1384;22.127839,203;22.127839,1378;53.981118,1384;40.559102,1378;23.242759,1384;21.682017,199;53.814322,1378;40.610481,1382;21.440680,1382;20.913617,1384;53.888045,1382;40.711450,1378;24.450980,203;69.088386,48;61.484103,1380;18.787536,201;18.041200,1384;18.304489,1382;19.617278,1378;18.304489,1380;54.121711,1384;24.061800,1384;21.185139,1382;21.051500,1382;21.051500,1378;20.913617,1378;54.112197,205;21.563025,1384;19.617278,1384;21.185139,1380;21.314789,1380;36.799045,1384;53.932594,1384;23.075702,1384;38.591158,199;21.185139,1380;19.424227,1380;21.185139,1380;53.940068,1378;22.989700,1378;33.944880,1378;38.869054,1384;22.334685,1382;21.797836,1380;53.812023,1378;18.787536,199;69.658583,49;61.051038,1384;19.979400,1380;19.010300,205;19.424227,1380;19.802112,1384;19.802112,1384;21.185139,1382;54.100173,1384;20.149733,1382;21.440680,1384;22.334685,205;21.051500,1380;21.051500,1382;21.051500,1378;54.045279,1382;20.913617,1378;36.195317,1378;21.563025,205;21.051500,1382;35.661417,1378;21.185139,1382;54.107432,1378;21.563025,1382;22.334685,1378;21.051500,203;22.232493,1382;21.051500,1378;32.085260,1380;54.049568,1384;22.812412,1384;21.185139,1378;21.563025,199;27.038037,1378;20.771213,1382;20.771213,1380;54.161750,1382;22.532125,1382;20.913617,1384;21.440680,201;21.185139,1380;21.051500,1378;21.314789,1382;54.056843,1384;29.838154,1378;20.913617,1378;20.771213,199;20.771213,1384;21.314789,1378;36.297895,1380;53.969141,1384;23.634280,1382;22.434527,1380;21.314789,203;69.496998,48;70.764129,50;70.781565,47;70.765583,47;70.764219,50;70.770015,49;70.782246,43;70.769062,45;70.769802,48;70.774975,50;70.766909,49;70.768461,52;70.759449,49;69.032628,24;83.075887,104];
X=Z(1:132,end-1);
Y=Z(1:132,end);

% Train up GMM on this data
chmm=struct('hmmone',[],'hmmtwo',[]);
chmm.hmmone.K=2; 
chmm.hmmtwo.K=2;

%chmm=chmmrandinit(X,Y,chmm,['full';'full']);
chmm=chmminit(X,Y,chmm,['full';'full']);
% $$$ chmm.hmmone.state(1).Mu=simchmm.hmmone.state(2).Mu';
% $$$ chmm.hmmone.state(1).Cov=simchmm.hmmone.state(2).Cov;
% $$$ chmm.hmmtwo.state(1).Mu=simchmm.hmmtwo.state(1).Mu';
% $$$ chmm.hmmtwo.state(1).Cov=simchmm.hmmtwo.state(1).Cov;
% $$$ chmm.hmmone.state(2).Mu=simchmm.hmmone.state(1).Mu';
% $$$ chmm.hmmone.state(2).Cov=simchmm.hmmone.state(1).Cov;
% $$$ chmm.hmmtwo.state(2).Mu=simchmm.hmmtwo.state(2).Mu';
% $$$ chmm.hmmtwo.state(2).Cov=simchmm.hmmtwo.state(2).Cov;


disp('Means of CHMM initialisation');
chmm.hmmone.state(:).Mu;
chmm.hmmtwo.state(:).Mu;

% Train up HMM on observation sequence data using Baum-Welch
% This uses the forward-backward method as a sub-routine
disp('We will now train the CHMM using Baum/Welch');
disp(' ');
disp('Press a key to continue');
pause
disp('Estimating CHMM');

chmm.train.cyc=100;
chmm.hmmone.obsmodel='Gauss';
chmm.hmmtwo.obsmodel='Gauss';
chmm.hmmone.train.obsupdate=ones(1,chmm.hmmone.K);% update observation models ?
chmm.hmmtwo.train.obsupdate=ones(1,chmm.hmmtwo.K);% update observation models ?
chmm.hmmone.train.init=1;         % Yes, we've already done initialisation
chmm.hmmtwo.train.init=1;         % Yes, we've already done initialisation


chmm=chmmtrain(X,Y,T,chmm);
disp('Means');
%chmm.hmmone.state.Mu
chmm.hmmtwo.state.Mu


disp('Initial State Probabilities, Pi');
chmm.Pi
disp('State Transition Matrix, P');
chmm.P

% For lag-1 models,  decoding is done identically to 
% decoding in standard HMMs with larger state space
disp('We will estimate the Viterbi path of the CHMM, ');
disp(' ');
disp('Press a key to continue');
pause

chmm.P=chmm.P';                        % Will's code uses transpose 
% joint viterbi path
[block.joint]=chmmdecode(X,Y,T,chmm);

% We are often interested in the Viterbi path of the individual chains. We can 
% estimate these in two ways: a) marginalise the transition probabilities and 
% run decoding for each of the two HMMs; b) marginalise the joint state
% sequence directly 
disp('We will estimate the Viterbi path of the each HMM, ');
disp(' ');
disp('First obtain Viterbi path directly from joint state sequence');
disp('Press a key to continue');
pause


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
  



disp(['Second obtain Viterbi path by decoding each of the HMMs using a' ...
      ' marginalised transition propability']);
disp('Press a key to continue');
pause

% single chain viterbi path
  % priors for individual chains
  Pone=ones(1:chmm.hmmone.K);
  Pone=Pone./length(Pone);
  Ptwo=ones(1:chmm.hmmtwo.K);
  Ptwo=Ptwo./length(Ptwo);

  % first chain; integrate out 2nd chain
  hmm=chmm.hmmone;
  hmm.P=mdsum(mdprod(hmm.P,Ptwo,3),3)';
  [block.hmmone]=hmmdecode(X,T,hmm);

  % second chain; integrate out 1st chain
  hmm=chmm.hmmtwo;
  hmm.P=mdsum(mdprod(hmm.P,Pone,2),2)';
  [block.hmmtwo]=hmmdecode(Y,T,hmm);

  % find incorrect classification
  inccla=find((diff(data.jointclass(1:T))~=0)-(diff(block.joint.q_star)~=0))+1;
  Tvec=1:T;
  
  figure;
  subplot(611),plot(Tvec,data.jointclass(1:T)),title('Input Viterbi');
  axis off
  subplot(613),plot(Tvec,block.joint.q_star),
  title('Joint Viterbi');
  axis off

