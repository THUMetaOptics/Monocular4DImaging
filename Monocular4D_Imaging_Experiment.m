%% Image preprocessing
clear,close all
I1 = imread('raw.bmp');

cx=-100;
cy=-100;
ax=-226;
ay=23;

result1 = I1(698+cx:1698+cx,2087+cy:3987+cy);
result2 = I1(3525+cx+ay:4525+cx+ay,2298+ax+cy:4198+ax+cy);

ra = 0.46;
xc = 849*ra*0.53;
yc = 1800*ra*0.45;
 result1 = imresize(result1,ra,'nearest');
 result2 = imresize(result2,ra,'nearest');
sz = size(result1);
res1 = medfilt2(result1,[3,4]);
res2 = medfilt2(result2,[3,4]);

result1 = imadjust(res1);
result2 = imadjust(res2);
result1 = imadjust(result1);
result2 = imadjust(result2);
        
figure(1),subplot(1,2,1);imagesc(result1);colormap('gray');title('PSF1(x-polarized)')
subplot(1,2,2);imagesc(result2);colormap('gray');title('PSF2(y-polarized)')

%% Object segmentation and labeling
tic

BW2 = edge(result2,'Canny',0.08,1.6);%
BW1 = edge(result1,'Canny',0.08,1.6);

 se2 = [0 1 1 1 0; 1 1 1 1 1; 1 1 1 1 1; 1 1 1 1 1; 0 1 1 1 0];
se=strel('disk',1);
se2=strel('square',3);
 BW1=imdilate(BW1,se2);
 BW2=imdilate(BW2,se2);
BW1x=~BW1;
BW2x=~BW2;

figure(2),subplot(1,2,1);imagesc(BW1x);colormap('gray');title('Canny operator edge detection1')
subplot(1,2,2);imagesc(BW2x);colormap('gray');title('Canny operator edge detection2')

L01=bwlabel(BW1x,4);
L02=bwlabel(BW2x,4);
BI1=imbinarize(result1,'adaptive');
BI2=imbinarize(result2,'adaptive');
figure(3),subplot(1,2,1);imagesc(BI1);title('Label1')
subplot(1,2,2);imagesc(BI2);title('Label2')

IT1 = imbinarize(L01,1.7);
IT2 = imbinarize(L02,1.7);

L1=bwlabel(IT1,8);
L2=bwlabel(IT2,8);

figure(4),subplot(1,2,1);imagesc(L1);
subplot(1,2,2);imagesc(L2);

  IT1=imdilate(L1,se2);
 IT2=imdilate(L2,se2);
  ds=round(30000*ra^2);
 IT1 = bwareaopen(IT1,ds);
 IT2 = bwareaopen(IT2,ds);

 IT1=~IT1;
 IT1 = bwareaopen(IT1,ds);
 IT1=~IT1;
 IT2=~IT2;
 IT2 = bwareaopen(IT2,ds);
 IT2=~IT2;
 L1=bwlabel(IT1,8);
L2=bwlabel(IT2,8);
figure(5),subplot(1,2,1);imagesc(L1);
subplot(1,2,2);imagesc(L2);

%% SAD matching and depth map
stats1 = regionprops(L1,'Centroid','MajorAxisLength');
cens1 = cat(1, stats1.Centroid);
le1=cat(1, stats1.MajorAxisLength);

nl=size(le1);

depth=zeros(sz(1),sz(2));
I01=zeros(sz(1),sz(2));
I02=I01;
for i = 1:nl
         cx=cens1(i,2);
         cy=cens1(i,1);
         dx=-((cy/sz(2)-0.55)+abs(cy/sz(2)-0.5))^1*1*ra-abs((cy/sz(2)-0.50)-abs(cy/sz(2)-0.50))^2.3*8*ra;  
         dy=(0.5-(cy/sz(2))+abs(cy/sz(2)-0.5))^1*9*ra+(-(0.5-(cy/sz(2))-abs(cy/sz(2)-0.5)))^2.3*3*ra;          %Correction of optical imaging distoration

         cx=round(cx);
         cy=round(cy);
        sr=round((le1(i)+5)/2);
        SL1=L1;
        SL2=L2;
        SL1(SL1<i)=0;
        SL1(SL1>i)=0;
        SL2(SL2<i)=0;
        SL2(SL2>i)=0;
        obj_temp1 =  SL1(cx-sr+1:cx+sr,cy-sr+1:cy+sr);
        er=zeros(24);
        
        for ii = 1:30                                  
            for jj = 1:30
                if sqrt((ii-15)^2+(jj-15)^2)>=2&&sqrt((ii-15)^2+(jj-15)^2)<=9.3 
                    obj_temp2 =  SL2(cx-sr+1+ii-15:cx+sr+ii-15,cy-sr+1+jj-15:cy+sr+jj-15);
                    
                    er(ii,jj) = sum(sum(abs(obj_temp1-obj_temp2)));      
                              
                    
                else er(ii,jj)=9e19;
                    
                end
            end
        end
        era{i}=er;
        abs_min=min(min(er));                                  %SAD matching
        [y,x]=find(er==abs_min);                             
        x=x(1);
        y=y(1);
        if er(y,x-1)<9e10&&er(y,x+1)<9e10
            px=((er(y,x-1)-er(y,x))-(er(y,x+1)-er(y,x)))/((er(y,x-1)-er(y,x))+(er(y,x+1)-er(y,x)));
        else
            px=0;
        end
        if er(y-1,x)<9e10&&er(y+1,x)<9e10
            py=((er(y-1,x)-er(y,x))-(er(y+1,x)-er(y,x)))/((er(y-1,x)-er(y,x))+(er(y+1,x)-er(y,x)));
        else
            py=0;
        end

        x=x-15;
        y=y-15;
      
        x=x+0.5*px+dy;            
        y=y+0.5*py+dx;                                           % Sub pixel refinement
        yy=round(x/2);
        xx=round(y/2);
        cer{i}=er;
        
        if x==0
            if y>0
                angle=90;           
            else angle=270;               
            end           
        else
            if y>=0
                if x>0
                    angle=atand(y/x);
                else
                    angle=180-atand(-y/x);
                end                
            else if x>0
                    angle=360-atand(-y/x);
                else angle=180+atand(y/x);
                end
            end
        end                                                             %Convert coordinates to angles
        
        angle=360-angle;                      
        
        ag(i)=angle;
        zf0=212;
        id=1/(1/20-1/zf0);
        nd=sqrt((cy-yc)^2+(cx-xc)^2)/ra;
        idc=sqrt((3.88*nd/1000)^2+id^2);                %Correction of PSF for off-axis objects
        zf=20*idc/(idc-20);                 
        z(i)=1/((28500+zf*(angle-180))/28500*1000/zf);
        
        [r1,c1] = find(L1==i);
        [r2,c2] = find(L2==i);
        for j=1:length(r1)
        
        I01(r1(j)+xx,c1(j)+yy) = res1(r1(j),c1(j));
        end
        for j=1:length(r2)
          depth(r2(j)+xx,c2(j)+yy)=z(i);  
        I02(r2(j)-xx,c2(j)-yy) = res2(r2(j),c2(j));
        end
end

toc 
depth(depth==0)=NaN;
depth=depth(:,110:770)*100;

 figure(6),imagesc(depth);title('Depth map');set(imshow(depth),'alphadata',~isnan(depth)*0.95+0.05);set(gca,'linewidth',0.8);colormap(brewermap2(168,'*BuGn'));
 colorbar;caxis([0.25,0.30]*100);title('z');

 %% 2D Intensity image and polarization contrast
 I11=medfilt2(I01(:,110:770),[5,5]);
 I22=medfilt2(I02(:,110:770),[5,5]);
 figure(7),subplot(1,2,1);imshow(I11,[0,255]);colormap('gray');
subplot(1,2,2);imshow(I22,[0,255]);colormap('gray');

pr=imdivide(I22,I11);
sz2=size(pr);

for i=1:sz2(1)
    for j=1:sz2(2)
        if pr(i,j)==Inf
            pr(i,j)=NaN;
        end
    end
end
pr(pr==0)=NaN;
pr=medfilt2(pr,[5,5]);
maxp=max(max(pr));
minp=min(min(pr));
figure(8),imshow(pr,[]);title('x-p/y-p');set(imshow(pr),'alphadata',~isnan(pr)*0.92+0.08);colormap(brewermap2(168,'*BuPu'));caxis([minp,maxp]);