function xPhys= mainMulti(nelx,nely,volfrac,initialDesign,posttreat,problem)
%used to obtain compliance and time of our method

%volfrac:   global volume fraction
%nelx:   number of cells in horizontal direction
%nely:   number of cells in vertical direction
%initialDesign:   initial guess : "volfrac" for uniform density, "top88" for mono-scale optimization
%posttreat:   0: no post-treatment (faster); 1: posttreatment(better design)
%problem:   'MBB', 'Lshape' or 'Canti'

volfraclim=0.09;
volfraclim1=0.599;
volfraclim2=0.8;

if volfrac < volfraclim
macroTic=tic;
[cTheor,xdens,xcos,xsin,xcub]=topMulti(nelx,nely,volfraclim,initialDesign,problem); %macro-optimization
macroTime=toc(macroTic);
totdesignTic=tic;
[xPhys]=totalDesign(xdens,xcos,xsin,xcub,nelx,nely,volfrac,posttreat,problem); %get total design %à laisser le volufraction et pas le volumfrac_limit
totdesignTime=toc(totdesignTic);

% %print final design
figure(5)
colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;

elseif volfrac > volfraclim1 && volfrac < volfraclim2
    
macroTic=tic;
[cTheor,xdens,xcos,xsin,xcub]=topMulti(nelx,nely,volfraclim1,initialDesign,problem); %macro-optimization
macroTime=toc(macroTic);
totdesignTic=tic;
[xPhys]=totalDesign(xdens,xcos,xsin,xcub,nelx,nely,volfrac,posttreat,problem); %get total design
totdesignTime=toc(totdesignTic);

% % print final design
figure(5)
colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;

else
macroTic=tic;
[cTheor,xdens,xcos,xsin,xcub]=topMulti(nelx,nely,volfrac,initialDesign,problem); %macro-optimization
macroTime=toc(macroTic);
totdesignTic=tic;
[xPhys]=totalDesign(xdens,xcos,xsin,xcub,nelx,nely,volfrac,posttreat,problem); %get total design
totdesignTime=toc(totdesignTic);

% % print final design
figure(5)
colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;

end

% % % % print final design
% figure(5)
% colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;

%evaluate final design
compliance=evaluateTotalDesign(xPhys,3,problem)
totalTime=macroTime+totdesignTime






