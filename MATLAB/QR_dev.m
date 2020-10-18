clc;clear;close all;


%% Setup

A=binread('A.bin');


%%

[U,S,V]=svd(A,'econ');


%% 
U_pdp=binread('U.bin');
P=binreadint('Pall.bin');

%% Check U
plot(U(:,1))
hold on
plot(U_pdp(:,1),'ro')
%%
close
%%
P=[P,P+length(U_pdp)/2];

%%
check=3

Usamp=U((P+1),:);
Usamp_pdp=binread('Usamp.bin');

Usamp_err= abs(Usamp)-abs(Usamp_pdp);

subplot(2,1,1)
plot(Usamp(:,check))
hold on
plot(Usamp_pdp(:,check),'r')
subplot(2,1,2)
plot(Usamp_err(:,check))
%%
close
%%
[UU,SS,VV]=svd(Usamp,'econ');

pinvcheck=VV'*diag(1./diag(SS))*UU';

%%
pinvUsamp=pinv(Usamp);

plot(pinvUsamp(check,:))
hold on
plot(pinvcheck(check,:),'r')
%%
pinvUsamp_pdp=binread('pinvUsamp.bin');

pinvUsamp_err= abs(pinvUsamp)-abs(pinvUsamp_pdp);

subplot(2,1,1)
plot(pinvUsamp(check,:))
hold on
plot(pinvUsamp_pdp(check,:),'r')
subplot(2,1,2)
plot(pinvUsamp_err(check,:))








%%
numModes=100;

Us=U(:,1:numModes);

[~,~,P]=qr(Us',0);
P=P(1:numModes)

%%
nPoints=20;
assert(nPoints>numModes)
assert(nPoints< size(A,1))


%% ODIEM E

while length(P)<nPoints
    lP=length(P);
    Usamp=zeros(lP,size(Us,2));
    for i=1:lP
        Usamp(i,:)=Us(P(i),:);
    end
    [~,~,V]=svd(Usamp,'econ');
    VT=V';
    
    r= (VT(end,:)*Us').^2;
    [~,I]=sort(r)
    I=flip(I)
    
    % add to
    for i=1:length(I)
        found=0;
        for j=1:length(P)
            if(P(j)==I(i))
                found=1;
                break
            end
        end
        if(found==0)
            P=[P,I(i)]
            break;
        end
    end
end

%%
lP=length(P);
Usamp=zeros(lP,size(Us,2));
for i=1:lP
    Usamp(i,:)=Us(P(i),:);
end








