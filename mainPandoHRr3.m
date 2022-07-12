% parpool('local',10)
density=0.01:0.001:1;
m=size(density,2);

% 
% parfor j=1:m
% for j =1:10
for j =800:900
    locfilename=['pointHRr3',num2str(j)];
    if ~isfile(locfilename)
        [v0,c1,v1,c2,v2]=mainTop88r3(20,10,density(j));
        locfileID=fopen(locfilename,'w');
        fprintf(locfileID,'%14.10f %14.10f %14.10f %14.10f %14.10f',v0,c1,v1,c2,v2);
    end
end
    
fclose('all');

% parfor j =1:m
% [v0,c1,v1,c2,v2]=mainTop88r3(20,10,density(j));
%     locfilename=['pointHr3',num2str(j)];
%     locfileID=fopen(locfilename,'w');
% fprintf(locfileID,'%14.10f %14.10f %14.10f %14.10f %14.10f',v0,c1,v1,c2,v2);
% end
%     
% fclose('all');