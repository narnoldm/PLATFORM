clc;clear all;close all;

%%


%Load in single snapshot to check format and also for use as centering
%subtraction


A1=binread('./POD/inj_bot.podqcons00000251');

%get vector length
l=length(A1)

%calc number of points
vars=10
cells=l/vars


fileStart=251
fileEnd=750
fileSkip=1
filePrefix='./POD/inj_bot.podqcons'


%quick plot check on variables
figure('Position',[10,10,900,1500])

for i=1:vars
    i;
    
    subplot(vars,1,i);
    plot(A1(((i-1)*cells+1):((i)*cells)))
    set(gca,'FontSize',10)
    set(gca,'LineWidth',2)
    ylabel(strcat('Variable',{' '},int2str(i)));
    
end
xlabel('Point ID')
saveas(gca,'Raw_Primitive_Variables.png')


%% Serial Scan over fileset to calculate normalization values 
% Step 2 in PDF





fileCount=0;
for fileId=fileStart:fileSkip:fileEnd
    fileCount=fileCount+1;
end
fileCount

varSums=zeros(vars,1);

% L2 Norm type
for fileId=fileStart:fileSkip:fileEnd
    fileId
    Asnap=binread(strcat(filePrefix,num2str(fileId,'%08d')));
    % Using snapshot 1 centering
    Apert=Asnap-A1;
    for i=1:vars
        if(i==2)
            varSums(i)= varSums(i)+  sum(Apert(((i-1)*cells+1):((i)*cells)).^2 +  Apert(((i)*cells+1):((i+1)*cells)).^2 +  Apert(((i+1)*cells+1):((i+2)*cells)).^2);
        else
            varSums(i)= varSums(i)+ sum(Apert(((i-1)*cells+1):((i)*cells)).^2);
        end
    end
end
varSums(3)=varSums(2)
varSums(4)=varSums(2);
varNorms=varSums./(cells*fileCount);
varNorms=sqrt(varNorms)














