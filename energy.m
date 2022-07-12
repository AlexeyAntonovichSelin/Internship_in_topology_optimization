% modulus_young_min = xlsread('database.csv', 1, 'BJ31:BJ121');
% 
% [ndata, text, alldata] = xlsread('database.csv');

clear all; 
c=readtable('database.csv');
cnum=c(15:105,61:62);
cnum=table2array(cnum);
% cnum=str2double(cnum)
M = cnum;
% depl=-cnum(1:end,8)
% deplM=-cnum(1:end,27) %machine displacement
% force=-cnum(1:end,28)
density_max = M(:, 1);
modulus_young_min = M(:, 2);
% 
% density_max=str2double(density_max)
% modulus_young_min=str2double(modulus_young_min)

% figure (1);
% loglog(density_max, modulus_young_min, 'ko');

allmatsE = modulus_young_min.';
allmatsRo = density_max.';

allmatsRo_min = min(allmatsRo);
allmatsE_min = min(allmatsE);
allmatsRo_max = max(allmatsRo);
allmatsE_max = max(allmatsE);

paretoE = [];
paretoRo = [];

[Ro_min, indmin] = min(allmatsRo);
[E_min, intmin] = min(allmatsE);

paretoE = allmatsE;
paretoRo = allmatsRo;
[lol, indice] = sortrows(paretoRo');
paretoRo = paretoRo (indice);
paretoE = paretoE (indice);
%paretoRo = paretoRo (1:4);
%paretoE = paretoE (1:4);

i=2;
maxi=size(paretoRo,2);
while i<=maxi
    if paretoE(i) < paretoE(i-1)
        paretoE=[paretoE(1,1:i-1) paretoE(1,i+1:end)];
        paretoRo=[paretoRo(1,1:i-1) paretoRo(1,i+1:end)];
        i=i-1;
        maxi=maxi-1;
    end
    i=i+1;
end

i=2
maxim=size(paretoRo,2);
while i<maxim

    if paretoE(i)>paretoE(i-1) && paretoRo(i-1)==paretoRo(i)
        paretoE(i-1)=[];
        paretoRo(i-1)=[];
        i=i-1;
        maxim=maxim-1;
    elseif paretoE(i-1)>paretoE(i) && paretoRo(i-1)==paretoRo(i)
        paretoE(i)=[];
        paretoRo(i)=[];
        i=i-1;
        maxim=maxim-1;
        elseif paretoRo(i)>paretoRo(i-1) && paretoE(i-1)==paretoE(i)
        paretoE(i)=[];
        paretoRo(i)=[];
        i=i-1;
        maxim=maxim-1;
    end
    i=i+1;
end


i=size(paretoRo,2);
while i>1
    if paretoRo(i)>paretoRo(i-1) && paretoE(i-1)==paretoE(i)
        paretoE(i)=[];
        paretoRo(i)=[];
    elseif paretoE(i)>paretoE(i-1) && paretoRo(i-1)==paretoRo(i)
        paretoE(i)=[];
        paretoRo(i)=[];
    end
    i=i-1;
end

% figure (11);
% plot(density_max, modulus_young_min, 'r*')
% hold on
% plot(paretoRo, paretoE, 'ko')
% hold off

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

density=0.01:0.001:1;
density=density';
m=size(density,1);
% complfilename=['complHRr3.txt'];
complfilename=['1000complHRr3.txt'];
complfileID=fopen(complfilename);
compliance=textscan(complfileID,'%24.10f');
compliance=reshape(compliance{1,1},5,991)';
v0=compliance(:,1);
c1=compliance(:,2);
v1=compliance(:,3);
c2=compliance(:,4);
v2=compliance(:,5);
fclose('all')

cP1dec=c1;
for i=2:m
    cP1dec(i)= min(c1(i),cP1dec(i-1));
end

dif=max(c1-cP1dec);

fclose('all')
% 
figure(1)
plot(density,c1-cP1dec)
figure(2)
plot(density,c1)


%1st try : same value in d=0.5 and 1


% legend('raw','0.2','0.5','1','2','5')
% hold off

%%no longer every n but find through 1 point (or more?)

nlow=0.01;
nup=100;
% di=212 % density coordinate of intersection between data and model
di=100;
data=c1(di);
% while nup-nlow>0.01
%     n=(nup+nlow)/2;
%     model=1./density.*exp(-(density.^n)./n);
%     aaa=c1(end)/model(end);
%     modelAdjusted=aaa*model;
%     modeldi=modelAdjusted(di);
%     if modeldi > datadi
%         nup=n;
%     else
%         nlow=n;
%     end
% end
% n
while nup-nlow>0.01
    n=(nup+nlow)/2;
    model=n./density(2:end-1)+(data(end)-n).*density(2:end-1).^(n/(data(end)-n));
    modeldi=model(di);
    if modeldi > data
        nup=n;
    else
        nlow=n;
    end
end
n;

figure(12)
% plot(density,(modelAdjusted-c1)./c1)
plot(density(1:end-2),(model-c1(1:end-2))./c1(1:end-2))

figure(100)
% plot(density,modelAdjusted)
plot(density(1:end-2),model)
hold on
plot(density(1:end-2),c1(1:end-2))
% plot(density,c1)
hold off


% %%see influence of point position
% maxDif=[];
% allN=[];
% for di=1:1:100
%     nlow=0.01;
%     nup=100;
%     %di=4950 % density coordinate of intersection between data and model
%     datadi=c1(di);
%     while nup-nlow>0.01
%         n=(nup+nlow)/2;
%         model=1./density.*exp((density.^n)./n);
%         aaa=c1(end)/model(end);
%         modelAdjusted=aaa*model;
%         modeldi=modelAdjusted(di);
%         if modeldi > datadi
%             nup=n;
%         else
%             nlow=n;
%         end
%     end
%     maxdif=max(abs((modelAdjusted-c1)./c1));
%     maxDif=[maxDif, maxdif];
%     allN=[allN, n];
% end
% figure(13)
% plot(1:1:100,maxDif)
% figure(14)
% plot(1:1:100,allN)

F = 20000; %N
Umax = 0.005; %m
l = 2;%m
h = 0.5;%m
t=0.01;%m

% new_density=[];
% for i=1:size(modulus_young_min,1)
% %     new = exp(-(70368744177664*(modulus_young_min(i)*Umax*t/F)^(7036444963986473/70368744177664))/7036444963986473)/(modulus_young_min(i)*Umax*t/F)
% %       new = exp(-(70368744177664*(modulus_young_min(i))^(7036444963986473/70368744177664))/7036444963986473)/modulus_young_min(i)
% %       new = exp(-(2251799813685248*modulus_young_min(i))^(8941414099661619/2251799813685248))/8941414099661619)/modulus_young_min(i);
%       new = exp(-lambertw(0, 1/modulus_young_min(i)^n)/n)/modulus_young_min(i);
%       new_density=[new_density new];
% end
% 
% new_youngus=[];
% for i=1:size(modulus_young_min,1)
%       new_1 = 1/density_max(i)*exp((density_max(i)^n)/n);
%       new_youngus=[new_youngus new_1];
% end

% 
% figure(17)
% loglog(modulus_young_min, new_density, 'ko')
% figure(171)
% loglog(density_max, new_youngus, 'ko')
% % hold on-
% figure(18)
% plot(density_max, modulus_young_min, 'ko')
% % figure(19)
% plot(modulus_young_min, density_max, 'ko')
% figure(20)
% plot(density, c1, 'b-')
% hold off

% new_indice=[];
% for i=1:size(new_density,2)
%       new_2 = density_max(i)*exp(-lambertw(0, 1/(modulus_young_min(i)*Umax*t/F)^n)/n)/(modulus_young_min(i)*Umax*t/F);
%       new_indice=[new_indice new_2];
% end
% 
% figure(21)
% % loglog(new_density', new_indice, 'ko')
% plot(density_max, new_indice, 'ko')
% 
% new_indice1=[];
% for i=1:size(paretoRo,2)
%       new_3 = paretoRo(i)*exp(-lambertw(0, 1/(paretoE(i)*Umax*t/F)^n)/n)/(paretoE(i)*Umax*t/F);
%       new_indice1=[new_indice1 new_3];
% end
% 
% figure(22)
% % loglog(new_density', new_indice, 'ko')
% plot(paretoRo, new_indice1, 'ko')

% volfrac=0.01:0.01:1
% volfrac=volfrac'
% new_complisance=[];
% for i=1:size(volfrac, 1)
%       new_4 = 1/volfrac(i)*exp((volfrac(i)^n)/n);
%       new_complisance=[new_complisance new_4];
% end
% 
% 
% new_volfrac=[];
% for i=1:size(volfrac, 1)
% %       new_5 = exp(-lambertw(0, 1/(new_complisance(i)*Umax*t/F)^n)/n)/(new_complisance(i)*Umax*t/F);
%         new_5 = exp(-(2251799813685248*(new_complisance(i)*Umax*t/F)^(8941414099661619/2251799813685248))/8941414099661619)/(new_complisance(i)*Umax*t/F)
%         new_5 = exp(-lambertw(0, 1/(new_complisance(i)*Umax*t/F)^4)/4)/(new_complisance(i)*Umax*t/F);
%         new_volfrac=[new_volfrac new_5];
% end
% 
% figure (45)
% plot (volfrac, new_complisance, 'ko')
% hold on
% plot(new_volfrac, new_complisance, 'ro')
% hold off

f_1=[];
for i=1:size(modulus_young_min, 1)
% y=(1000000*modulus_young_min(84)*t*Umax)/F;
y=(1000000000*modulus_young_min(i)*t*Umax)/F;
xmin=0;
xmax=1;
while xmax-xmin> 1.0000e-06
   xmid=(xmax+xmin)/2;
   fmid=(1/xmid)*exp(-(xmid^n)/n)
   if fmid<y
      xmax=xmid;
   else
      xmin=xmid;
   end
end
xmid
f_1=[f_1 xmid]   
end

f_1Pareto=[];
for i=1:size(paretoE, 2)
y=(1000000000*paretoE(i)*t*Umax)/F;
xmin=0;
xmax=1;
while xmax-xmin> 1.0000e-06
   xmid=(xmax+xmin)/2;
   fmid=(1/xmid)*exp(-(xmid^n)/n)
   if fmid<y
      xmax=xmid;
   else
      xmin=xmid;
   end
end
xmid
f_1Pareto=[f_1Pareto xmid]   
end

figure (45)
plot (f_1',(1000000000*modulus_young_min.*t*Umax)/F, 'ko')
hold on
plot (f_1Pareto',(1000000000*paretoE.*t*Umax)/F, 'r+')
hold off

indice=[];
for i=1:size(modulus_young_min, 1)
okay = density_max(i)*f_1(i);
indice=[indice okay];
end


indicePareto=[];
for i=1:size(paretoE, 2)
okay = paretoRo(i)*f_1Pareto(i);
indicePareto=[indicePareto okay];
end

figure (46)
plot (density_max,indice, 'ko')
hold on
plot (paretoRo, indicePareto, 'r+')
hold off

figure (47)
loglog (density_max,indice, 'ko')
hold on
loglog (paretoRo, indicePareto, 'r+')
hold off

% volumefractionmat=[];
% for i=1:size(density_max,1)
% calcule=f_1(i)/density_max(i);
% volumefractionmat=[volumefractionmat calcule];
% end

% volumefractionmatPareto=[];
% % for i=1:size(paretoRo,2)
% for i=10:10
% % calculeP=f_1Pareto(i)/paretoRo(i);
% calculeP=f_1Pareto(10)/paretoRo(10);
% volumefractionmatPareto=[volumefractionmatPareto calculeP];
% end

figure (48)
plot (f_1,indice, 'ko')
hold on
plot (f_1Pareto, indicePareto, 'r+')
hold off


figure (49)
loglog (f_1,indice, 'ko')
hold on
loglog(f_1Pareto, indicePareto, 'r+')
hold off
