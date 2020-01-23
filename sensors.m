%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model estimation procedure %                                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load time series RR intervals
load RR
% load the Boolean time series of AF detections for RR
load D
% define window' size
w=30;
% compute number of windows of length w in the data set
lw=ceil(length(RR)/w);
pos_af=zeros(lw,1);
for k=1:lw
    % take a time series of RR interval of length w                
    data=RR(1+w*(k-1):k*w,1);
    % a window is defined to be in AF if number of AF (Atrial Fibrillation) detections in data is >=w/2
    % define pos_af as the boolean vector identifying the AF windows
    if length(find(D(1+w*(k-1):k*w)==1))>w/2
       pos_af(k)=1;
    end;
    
    % Computation of covariates
    
    % Descriptive metrics
        %mean
        meanRR=mean(data);
        %median
        medianRR=median(data);
        % Pearson coefficient of variation
        VRR=std(data)/mean(data)
        % Mean absolute dispersion with respect to the median
        VmeRR=sum(abs(data-medianRR))/medianRR;
    % Symbolic recurrence measures
    
        [A matrix SRR cv]=func_SRP(data,m);
        %SRR{i} is the symbolic recurrence rate to symbol i
        SRRT=[SRR{1} SRR{2} SRR{3} SRR{4} SRR{5} SRR{6}];
        
        % Compute the following symbolic recurrence quantification analysis
        % of symbolic recurrence plots of the increasing and decreasing
        % symbols as well as of the complete recurrence plot
  
        % RR:   Recurrence rate RR; The percentage of recurrence points in
        % an recurrence matrix (or recurrence plot)
        % DET:  Determinism DET; The percentage of recurrence points which form
        %       diagonal lines
        % ENTR: Entropy ENTR; The Shannon entropy of the probability distribution of the diagonal
        %       line lengths p(l)
        % L:    Averaged diagonal line length L; The average length of the diagonal lines
        % V: contains the distribution of vertical lines together with their entropy
        % V(1,3) Shannon entropy of the probability distribution of the
        % vertical line lengths p(l)
        % V(1,4) mean length of vertical lines
        
        % SRP of Increasing symbol
        [RRinc,DETinc,DEThatinc,ENTRinc,Linc,Lhatinc,Vinc,ddinc,ndinc]= Recu_SRQA(A{factorial(m)},2,0);
        % SRP of Decreasing symbol
        [RRdec,DETdec,DEThatdec,ENTRdec,Ldec,Lhatdec,Vdec,dddec,nddec] = Recu_SRQA(A{1},2,0);
        % SRP all symbols          
        [RRt,DETt,DEThatt,ENTRt,Lt,Lhatt,Vt,ddt,ndt] = Recu_SRQA(matriz,2,0);
        
        
        %Storage of all metrics in a matrix         
        measure(k,:)=[mean(data), median(data),std(data)/mean(data),...
                        sum(abs(data-median(data)))/median(data),...
                        lisfrec,DETt,ENTRt,Vinc(1,3),Vdec(1,3),Vinc(1,4),Vdec(1,4)];
end;

   %compute the logistic model   
    mod=fitglm(measure,pos_af,'Distribution','binomial','Link','logit');
   %compute the estimated probability of AF
    prob=mod.Fitted.Probability;
   
    
   % compute the threshold for AF classification
   Tlist=[0.001:0.001:1];
        dista=zeros(length(Tlist),1);
        for j=1:length(Tlist)
            ths=Tlist(j);
            % true positive
        TP=length(find(prob(find(pos_af>0))>ths));
            % true negative
        TN=length(find(prob(find(pos_af==0))<ths));
            % false negative
        FN=length(find(prob(find(pos_af>0))<ths));
            % false positive
        FP=length(find(prob(find(pos_af==0))>ths));
        %sensitivity
        TPR=TP/(TP+FN);
        FPR=FP/(FP+VN);
        %specificity
        SPC=1-FPR;
        xc(j)=FPR;
        yc(j)=TPR;
        dista(j)=FPR^2+(1-VPR)^2;
        end;
        
        % threshold
        ths=Tlist(min(find(dista==min(dista))));
   
        
 %%% Functions used in the code above
 
 function SRR =func_SRP(x,m)

 % x data set
 % m is the embbeding dimension
T=length(x);

for i=1:T-(m-1);
    [B,mh(i,:)]=sort(x(i:d:i+d*(m-1))');
end;    

symbs=perms(1:m);

for i=1:factorial(m)
    loc{i}=ismember(mh,symbs(i,:),'rows');
    A{i}=kron(loc{i},loc{i}');
    SRR{i}=sum(sum(A{i}))/(T-(m-1))^2;
end;


function [RR,DET,DEThat,ENTR,L,Lhat,V,dd,nd] = Recu_SRQA(RP,Lmin,I)
% Recurrence quantification analysis of symbolic recurrence plots
% RP:  the Symbolic Recurrence Plot
% I:   the indication marks (I=0 RP is the symmetry matrix
%                            I=1 RP is the asymmetry matrix)  
% RR:   Recurrence rate RR; The percentage of recurrence points in an RP
%       Corresponds to the correlation sum;
% DET:  Determinism DET; The percentage of recurrence points which form
%       diagonal lines
% ENTR: Entropy ENTR; The Shannon entropy of the probability distribution of the diagonal
%       line lengths p(l)
% L:    Averaged diagonal line length L; The average length of the diagonal lines
% V: contains the distribution of vertical lines together with their entropy
% see below
%% 
%Lmin=2;
if nargin < 2
    I=0;
end
N1=size(RP,1);
Yout=zeros(1,N1);
for k=2:N1
    On=1;
    while On<=N1+1-k
        if RP(On,k+On-1)==1
            A=1;off=0;
            while off==0 & On~=N1+1-k
                if RP(On+1,k+On)==1
                    A=A+1;On=On+1;
                else
                    off=1;
                end
            end
            Yout(A)=Yout(A)+1;
        end 
        On=On+1;
    end
end
if I==0
    S=2*Yout;
end       
if I==1
    RP=RP';
    for k=2:N1
        On=1;
        while On<=N1+1-k
            if RP(On,k+On-1)==1
                A=1;off=0;
                while off==0 & On~=N1+1-k
                    if RP(On+1,k+On)==1
                        A=A+1;On=On+1;
                    else
                        off=1;
                    end
                end
                Yout(A)=Yout(A)+1;
            end 
            On=On+1;
        end
    end
    S=Yout;
end
%% calculate the recurrence rate (RR)
SR=0;
for i=1:N1
    SR=SR+i*S(i);
end
RR=SR/(N1*(N1-1));
%% calculate the determinism (%DET)
if SR==0
    DET=0;
else
    DET=(SR-sum(S(1:Lmin-1)))/SR;
end
DEThat=(SR-sum(S(1:Lmin-1)))/N1;
%% calculate the ENTR = entropy (ENTR)
pp=S/sum(S);
entropy=0;
F=find(S(Lmin:end));
l=length(F);
if l==0
    ENTR=0;
else
    F=F+Lmin-1;
    ENTR=-sum(pp(F).*log(pp(F)));
end
dd=S(F);
nd=1+[1:length(S)];
%% calculate Averaged diagonal line length (L)
if sum(S(Lmin:end))==0
    L=0;
else
L=(SR-sum([1:Lmin-1].*S(1:Lmin-1)))/sum(S(Lmin:end));
end;
Lhat=L/N1;

%% calculate vertical line length distributions (V)
VV=[];

for i=1:N1
    RPvi=RP(:,i);
     cv=cellfun(@abs,strsplit(char(RPvi'),char(0)),'uniform',false);
     v=zeros(length(cv),1);
     for j=1:length(cv)
         v(j)=sum(cv{j});
     end;
     VV=[VV;v(find(v>1))];
 end;
 if length(VV)<1
     VV=0;
 end;
[av bv sv]=unique(VV);
distV=hist(sv,length(av));
pV=distV/sum(distV);
pV(pV==0)=1;
entV=-sum(pV.*log(pV));
V=[av distV' repmat(entV,length(av),1) repmat((sum(av.*distV')/sum(distV)),length(av),1)];
