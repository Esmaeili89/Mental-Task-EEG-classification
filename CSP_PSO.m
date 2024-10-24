clc
cd 'E:\CODE&&dataset\--main code';
load eegdata;

class=zeros(1,25);
for i=1:13
i1=(i-1)*25;
    class(1,i1+1:i1+5)=1;
    class(1,i1+6:i1+10)=2;
    class(1,i1+11:i1+15)=3;
    class(1,i1+16:i1+20)=4;
    class(1,i1+21:i1+25)=5;
end

class_train=zeros(1,260);
class_train(1,1:52)=1;
class_train(1,53:104)=2;
class_train(1,105:156)=3;
class_train(1,157:208)=4;
class_train(1,209:260)=5;

class_test=zeros(1,65);
class_test(1,1:13)=1;
class_test(1,14:26)=2;
class_test(1,27:39)=3;
class_test(1,40:52)=4;
class_test(1,53:65)=5;

%%--------make window-------------------------------------------------------
wave=zeros(7,125,20,325);
for i=1:325
taske1=data{1,i}{1,4};
    for j=1:20   
       j1=125*j;        
    wave(:,:,j,i)=taske1(:,j1-124:j1);
    end
end

TP_win=zeros(1,20);
for win=1:20
R1=zeros(7,7);
k=0;
for p=1:13
    cc=(p-1)*25;
    for t=cc+1:cc+5
%         E1=data{1,t}{1,4}';
% wave1=wave(:,:,win,t);
wave1=wave(:,:,win,t);
          E1=wave1';
        sum=0;
        for y=1:125
        sum=sum+E1(y,:);   
        end
        sum=sum/125;
        Rt=zeros(7,7);
        for y=1:125
            r1=E1(y,:)-sum;
            r=r1'*r1;
            Rt=Rt+r;
        end
        cov=Rt/trace(Rt);
        R1=R1+cov;
        k=k+1;
    end
end
R1=R1/k;
%---------------------------**multiply**-----------------------------------
R2=zeros(7,7);
k=0;
for p=1:13
    cc=(p-1)*25;
    for t=cc+6:cc+10
        wave1=wave(:,:,win,t);
        E2=wave1';

        sum=0;
        for y=1:125
        sum=sum+E2(y,:);   
        end
        sum=sum/125;
        Rt=zeros(7,7);
        for y=1:125
            r2=E2(y,:)-sum;
            r=r2'*r2;
            Rt=Rt+r;
        end
        cov=Rt/trace(Rt);
        R2=R2+cov;
        k=k+1;
    end
end
R2=R2/k;
%---------------------^^letter^^-------------------------------------------
save_cov=zeros(65,7,7);
R3=zeros(7,7);
k=0;
s=zeros(65,125,7);
for p=1:13
    cc=(p-1)*25;
    for t=cc+11:cc+15
        wave1=wave(:,:,win,t);
        E3=wave1';

        sum=0;
        for y=1:125
        sum=sum+E3(y,:);  
        s(k+1,y,:)=sum;
        end
        sum=sum/125;
        Rt=zeros(7,7);
        for y=1:125
            r3=E3(y,:)-sum;
            r=r3'*r3;
            Rt=Rt+r;
        end
        cov=Rt/trace(Rt);
        save_cov(k+1,:,:)=cov;
          if k==34
              save_cov(k+1,:,:)=0;
          else
               R3=R3+cov;
%              nans = find(isnan(save_cov));
%              if ~isempty(nans)
%              save_cov(nans,:,:) = [];
%              end
        
             
          end

        
        SAVE_R3(k+1,:,:)=R3;
        k=k+1;
    end
end
R3=R3/k;
%---------------------^^Rotation^^-----------------------------------------

R4=zeros(7,7);
k=0;
for p=1:13
    cc=(p-1)*25;
    for t=cc+16:cc+20
        wave1=wave(:,:,win,t);
        E4=wave1';

        sum=0;
        for y=1:125
        sum=sum+E4(y,:);   
        end
        sum=sum/125;
        Rt=zeros(7,7);
        for y=1:125
            r4=E4(y,:)-sum;
            r=r4'*r4;
            Rt=Rt+r;
        end
        cov=Rt/trace(Rt);
        R4=R4+cov;
        k=k+1;
    end
end
R4=R4/k;

%2222222222222222222222222---counting---22222222222222222222222222222222222
R5=zeros(7,7);
k=0;
for p=1:13
    cc=(p-1)*25;
    for t=cc+21:cc+25
        wave1=wave(:,:,win,t);
        E5=wave1';

        sum=0;
        for y=1:125
        sum=sum+E5(y,:);   
        end
        sum=sum/125;
        Rt=zeros(7,7);
        for y=1:125
            r5=E5(y,:)-sum;
            r=r5'*r5;
            Rt=Rt+r;
        end
        cov=Rt/trace(Rt);
        R5=R5+cov;
        k=k+1;
    end
end
R5=R5/k;
%##########################################################################
total_R=R1+R2+R3+R4+R5;

[PC,V]=eig(total_R);
V=diag(V);
[junk,rindices]=sort(-1*V);
V=V(rindices);
PC=PC(:,rindices);
V1=zeros(size(V,1),size(V,1));
for g=1:size(V,1)
    V1(g,g)=V(g,1);
end
W=V1^(-0.5) *PC';
%@@@@@@@@@@@@@@@@@@@@@@@@--baseline--@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

Ra1=R2+R3+R4+R5;
Sa=W*R1*W';
Sa1=W*Ra1*W';
[PCa,Va]=eig((Sa+Sa')/2);
[PCa1,Va1]=eig((Sa1+Sa1')/2);

%--sorting eigVector && eigValue
Va=diag(Va);
[junk,rindices]=sort(-1*Va);
Va=Va(rindices);
PCa=PCa(:,rindices);
Vaa=zeros(size(Va,1),size(Va,1));
for g=1:size(Va,1)
    Vaa(g,g)=Vaa(g,1);
end


Va1=diag(Va1);
[junk,rindices]=sort(-1*Va1);
Va1=Va1(rindices);
PCa1=PCa1(:,rindices);
Vaa1=zeros(size(Va1,1),size(Va1,1));
for g=1:size(Va1,1)
    Vaa1(g,g)=Vaa1(g,1);
end


%--------------------------SFa=PCa*W;
Ua(:,1:2)=PCa(:,1:2);%m=2
Ua(:,3:4)=PCa(:,6:7);%m=2
% Pa=W'*Ua;
SFa=W'*Ua;
f1=zeros(65,4);

k=1;
for p=1:13
    cc=(p-1)*25;
    for t=cc+1:cc+5
        E1=data{1,t}{1,4}';
        Z1=SFa'*E1';
        sum=0;
        v1=zeros(1,4);
        for h=1:4
            v1(1,h)=var(Z1(h,:));
            sum=sum+v1(1,h);
        end
        for h=1:4
            f1(k,h)=v1(1,h)/sum;
            f1(k,h)=log(f1(k,h));
        end
   
        k=k+1;
    end
end
win_f1(:,(win-1)*4+1:win*4)=f1;
%$$$$$$$$$$$$$$$$$$$$--multipliction--$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
Rb1=R1+R3+R4+R5;
Sb=W*R2*W';
Sb1=W*Rb1*W';
[PCb,Vb]=eig((Sb+Sb')/2);
[PCb1,Vb1]=eig((Sb1+Sb1')/2);

%--sorting eigVector && eigValue
Vb=diag(Vb);
[junk,rindices]=sort(-1*Vb);
Vb=Vb(rindices);
PCb=PCb(:,rindices);
Vbb=zeros(size(Vb,1),size(Vb,1));
for g=1:size(Vb,1)
    Vbb(g,g)=Vbb(g,1);
end


Vb1=diag(Vb1);
[junk,rindices]=sort(-1*Vb1);
Vb1=Vb1(rindices);
PCb1=PCb1(:,rindices);
Vbb1=zeros(size(Vb1,1),size(Vb1,1));
for g=1:size(Vb1,1)
    Vaa1(g,g)=Vbb1(g,1);
end


%--------------------------SFa=PCa*W;
Ub(:,1:2)=PCb(:,1:2);%m=2
Ub(:,3:4)=PCb(:,6:7);%m=2
% Pa=W'*Ua;
SFb=W'*Ub;
f2=zeros(65,4);

k=1;
for p=1:13
    cc=(p-1)*25;
    for t=cc+6:cc+10
        E2=data{1,t}{1,4}';
        Z2=SFb'*E2';
        sum=0;
        v2=zeros(1,4);
        for h=1:4
            v2(1,h)=var(Z2(h,:));
            sum=sum+v2(1,h);
        end
        for h=1:4
            f2(k,h)=v2(1,h)/sum;
            f2(k,h)=log(f2(k,h));
        end
   
        k=k+1;
    end
end
win_f2(:,(win-1)*4+1:win*4)=f2;
%###########################--letter__#####################################
Rc1=R1+R2+R4+R5;
Sc=W*R3*W';
Sc1=W*Rc1*W';
[PCc,Vc]=eig((Sc+Sc')/2);
[PCc1,Vc1]=eig((Sc1+Sc1')/2);

%--sorting eigVector && eigValue
Vc=diag(Vc);
[junk,rindices]=sort(-1*Vc);
Vc=Vc(rindices);
PCc=PCc(:,rindices);
Vcc=zeros(size(Vc,1),size(Vc,1));
for g=1:size(Vc,1)
    Vcc(g,g)=Vcc(g,1);
end


Vc1=diag(Vc1);
[junk,rindices]=sort(-1*Vc1);
Vc1=Vc1(rindices);
PCc1=PCc1(:,rindices);
Vcc1=zeros(size(Vc1,1),size(Vc1,1));
for g=1:size(Vc1,1)
    Vcc1(g,g)=Vcc1(g,1);
end


Uc(:,1:2)=PCc(:,1:2);%m=2
Uc(:,3:4)=PCc(:,6:7);%m=2
% Pa=W'*Ua;
SFc=W'*Uc;
f3=zeros(65,4);

k=1;
for p=1:13
    cc=(p-1)*25;
    for t=cc+11:cc+15
        E3=data{1,t}{1,4}';
        Z3=SFc'*E3';
        sum=0;
        v3=zeros(1,4);
        for h=1:4
            v3(1,h)=var(Z3(h,:));
            sum=sum+v3(1,h);
        end
        for h=1:4
            f3(k,h)=v3(1,h)/sum;
            f3(k,h)=log(f3(k,h));
        end
   
        k=k+1;
    end
end
win_f3(:,(win-1)*4+1:win*4)=f3;
%###########################--Rotation--###################################
Rd1=R1+R2+R3+R5;
Sd=W*R4*W';
Sd1=W*Rd1*W';
[PCd,Vd]=eig((Sd+Sd')/2);
[PCd1,Vd1]=eig((Sd1+Sd1')/2);

%--sorting eigVector && eigValue
Vd=diag(Vd);
[junk,rindices]=sort(-1*Vd);
Vd=Vd(rindices);
PCd=PCd(:,rindices);
Vdd=zeros(size(Vd,1),size(Vd,1));
for g=1:size(Vd,1)
    Vdd(g,g)=Vdd(g,1);
end


Vd1=diag(Vd1);
[junk,rindices]=sort(-1*Vd1);
Vd1=Vd1(rindices);
PCd1=PCd1(:,rindices);
Vdd1=zeros(size(Vd1,1),size(Vd1,1));
for g=1:size(Vd1,1)
    Vdd1(g,g)=Vdd1(g,1);
end


Ud(:,1:2)=PCd(:,1:2);%m=2
Ud(:,3:4)=PCd(:,6:7);%m=2
% Pa=W'*Ua;
SFd=W'*Ud;
f4=zeros(65,4);

k=1;
for p=1:13
    cc=(p-1)*25;
    for t=cc+16:cc+20
        E4=data{1,t}{1,4}';
        Z4=SFd'*E4';
        sum=0;
        v4=zeros(1,4);
        for h=1:4
            v4(1,h)=var(Z4(h,:));
            sum=sum+v4(1,h);
        end
        for h=1:4
            f4(k,h)=v4(1,h)/sum;
            f4(k,h)=log(f4(k,h));
        end
   
        k=k+1;
    end
end
win_f4(:,(win-1)*4+1:win*4)=f4;
%###########################--Counting--###################################
Re1=R1+R2+R3+R4;
Se=W*R5*W';
Se1=W*Re1*W';
[PCe,Ve]=eig((Se+Se')/2);
[PCe1,Ve1]=eig((Se1+Se1')/2);

%--sorting eigVector && eigValue
Ve=diag(Ve);
[junk,rindices]=sort(-1*Ve);
Ve=Ve(rindices);
PCe=PCe(:,rindices);
Vee=zeros(size(Ve,1),size(Ve,1));
for g=1:size(Ve,1)
    Vee(g,g)=Vee(g,1);
end


Ve1=diag(Ve1);
[junk,rindices]=sort(-1*Ve1);
Ve1=Ve1(rindices);
PCe1=PCe1(:,rindices);
Vee1=zeros(size(Ve1,1),size(Ve1,1));
for g=1:size(Ve1,1)
    Vee1(g,g)=Vee1(g,1);
end


Ue(:,1:2)=PCe(:,1:2);%m=2
Ue(:,3:4)=PCe(:,6:7);%m=2
% Pa=W'*Ua;
SFe=W'*Ue;
f5=zeros(65,4);

k=1;
for p=1:13
    cc=(p-1)*25;
    for t=cc+21:cc+25
        E5=data{1,t}{1,4}';
        Z5=SFe'*E5';
        sum=0;
        v5=zeros(1,4);
        for h=1:4
            v5(1,h)=var(Z5(h,:));
            sum=sum+v5(1,h);
        end
        for h=1:4
            f5(k,h)=v5(1,h)/sum;
            f5(k,h)=log(f5(k,h));
        end
   
        k=k+1;
    end
end
win_f5(:,(win-1)*4+1:win*4)=f5;
end
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
test(1:13,:)=win_f1(1:13,:);
test(14:26,:)=win_f2(1:13,:);
test(27:39,:)=win_f3(1:13,:);
test(40:52,:)=win_f4(1:13,:);
test(53:65,:)=win_f5(1:13,:);
% test(26:65,:)=win_f5(26:65,:);

train(1:52,:)=win_f1(14:65,:);
train(53:104,:)=win_f2(14:65,:);
train(105:156,:)=win_f3(14:65,:);
train(157:208,:)=win_f4(14:65,:);
train(209:260,:)=win_f5(14:65,:);
% % % train(150,:)=0;%%-----------------56%
train(126,:)=0;%%-----------------56%

cd 'E:\COURSE';
load best_popu;
chrom1=zeros(110,20);%110 particle with 20 gene
for iu=1:110
    for ju=1:20
    chrom1(iu,ju)=round(rand(1));
    end
end

%------------number of windows that this chrom select----------------------
TP=zeros(30,110);%accuracy matrix with 110 particle and 30 epoch
num_particle=110;
x_gbest=struct('TP',0, 'num_chrom',1,'epoch',1, 'chromosome',[],'particle',[]);
x_pbest(3)=struct('TP',0, 'index',0,'particle',[]);
for u1=1:num_particle
    x_pbest(u1).TP=0;
end
%66666666666666666666666666666666666666666666666666666666666666666666666666


window_lenght=zeros(1,110);%show the number of 1 in each chrom and consider it as lenght 
index=zeros(110,20);

for i=1:110% i th particle    
    k=1;
    for j=1:20
        if chrom1(i,j)==1
            index(i,k)=j;
            k=k+1;
            window_lenght(1,i)=window_lenght(1,i)+1;
        end
    end
end
%444444444444444444444444444444444--PSO--4444444444444444444444444444444444
v1=zeros(110,4,325);%vi(t) for PSO
v2=zeros(110,4,325);%vi(t+1) for PSO
wmax=0.9;
wmin=0.4;
itmax=30; %Maximum iteration number
c1=1.4;
c2=1.4;
for epoch=1:itmax
w(epoch)=wmax-((wmax-wmin)/itmax)*epoch;
end
x_total_particle=zeros(110,4,325);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
for epoch=1:30

for num_chrom=1:110    %num_chrom: number of particles

%------------------------make training data--------------------------------
 ga_by_csp_data_train=zeros(260,4*window_lenght(1,num_chrom));
 for pc=1:window_lenght(1,num_chrom)
    I=index(num_chrom,pc);
    ga_by_csp_data_train(:,(pc-1)*4+1:pc*4)=train(:,(I-1)*4+1:I*4);
 end
 %------------------------make test data-------------------------------
 ga_by_csp_data_test=zeros(65,4*window_lenght(1,num_chrom));
 for p=1:window_lenght(1,num_chrom)
    I1=index(num_chrom,p);
    ga_by_csp_data_test(:,(p-1)*4+1:p*4)=test(:,(I1-1)*4+1:I1*4);
 end
 
 %----------------------*LDA*---------------------------------------

c = classify(ga_by_csp_data_test,ga_by_csp_data_train,class_train);

 tp=0;
 for t=1:65

     if(c(t,1)==class_test(1,t))
         tp=tp+1;
     end
 end
 TP(epoch,num_chrom)=tp;
 
 if x_pbest(num_chrom).TP<TP(epoch,num_chrom)
    x_pbest(num_chrom).TP=TP(epoch,num_chrom);
    x_pbest(num_chrom).index=epoch;
    x_pbest(num_chrom).particle=x_total_particle(num_chrom,:,:);
    
end%if    

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%44444444444444444444444444444444444444444444444444444444444444444444444444
end%...th particle
max_TP=TP(epoch,1);
max_index=1;
for u=1:110
    if TP(epoch,u)>max_TP
        max_TP=TP(epoch,u);
        max_index=u;
    end
end
        
if x_gbest.TP<max_TP
   x_gbest.TP=max_TP;
   x_gbest.num_chrom=max_index;
   x_gbest.epoch=epoch;
   x_gbest.chromosome=chrom1(max_index,:);
   x_gbest.particle=x_total_particle(max_index,:,:);
    
end%if    


for num_chrom=1:110

p1=x_pbest(num_chrom).particle;
g1=x_gbest.particle;
x1=x_total_particle(num_chrom,:,:);
v2(num_chrom,:,:)=w(epoch)*v1(num_chrom,:,:)+c1*rand(1)*(p1-x1)+c2*rand(1)*(g1-x1);
x_total_particle(num_chrom,:,:)=x_total_particle(num_chrom,:,:)+v2(num_chrom,:,:);

end
v1=v2;

end
%####################################################

for vv=1:num_chrom
    
  PSO_TP(vv)=x_pbest(vv).TP;
    
end
x_pbest_tp=PSO_TP/65*100;
x_gbest_tp=x_gbest.TP/65*100;

plot(x_gbest_tp);
xlabel('x_pbest');
xlabel('TP');

hold on;
plot(x_gbest_tp,'o');
hold on;
plot(x_gbest_tp,'+');

