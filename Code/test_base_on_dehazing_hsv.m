clc;
clear;
for z_p=1:12
I=imread(['images of project\set 1\blister_mixed_',num2str(z_p,'%02d'),'.png']);
% I=imread(['images of project\set 2\bl-130-',num2str(z_p,'%02d'),'.tiff']);

% I=imread('images of project\set 1\template.png');
% I=imread('images of project\set 2\template.tiff');
%I=imread(['images of project\set 1\blister_mixed_',num2str(11,'%02d'),'.png']);
% I=imread(['images of project\set 2\bl-130-',num2str(3,'%02d'),'.tiff']);

I_de=dehazing_main(I);

I_gray=im2gray(I);
se=ones(80,80);
I_close=imclose(I_gray,se);
I_open=imopen(I_close,se);
%     figure(z_p),subplot(221),imshow(I_gray),title('原始灰度图片');
%     subplot(222),imshow(I_open),title('闭开运算后的灰度图片');
I_bw=imbinarize(I_open,0.5);
%     subplot(223),imshow(I_bw),title('二值化后的灰度图片');
%统计标注连通域
%使用外接矩形框选连通域，并使用形心确定连通域位置
[l,b_n] = bwlabel(I_bw,8);
status_box_b=regionprops(l,'BoundingBox');
centroid_box = regionprops(l,'Centroid');
%     subplot(224),imshow(I_gray),title('板标记图片');hold on;
%     for i=1:b_n
%         rectangle('position',status_box(i).BoundingBox,'edgecolor','r');
%         text(centroid_box(i,1).Centroid(1,1)-15,centroid_box(i,1).Centroid(1,2)-15, num2str(i),'Color', 'r') 
%     end
%获取矩形中心点
center_xy=centroid_box(1).Centroid;

%获取矩形宽高
BoundingBox=status_box_b(1).BoundingBox(1,3:4)-[20,50];

%通过矩形中心坐标以及宽高获取每个分数的图像

I_Board=I_de(uint16(center_xy(1,2)-BoundingBox(1,2)/2):uint16(center_xy(1,2)+BoundingBox(1,2)/2),uint16(center_xy(1,1)-BoundingBox(1,1)/2):uint16(center_xy(1,1)+BoundingBox(1,1)/2),:);
I_Board_nor=I(uint16(center_xy(1,2)-BoundingBox(1,2)/2):uint16(center_xy(1,2)+BoundingBox(1,2)/2),uint16(center_xy(1,1)-BoundingBox(1,1)/2):uint16(center_xy(1,1)+BoundingBox(1,1)/2),:);

I_Board_hsv=rgb2hsv(I_Board);
I_Board_hsv(:,:,3)=I_Board_hsv(:,:,3)-I_Board_hsv(:,:,3)+0.5;
figure(1),subplot(221),imshow(I_Board_hsv),title('胶囊板图像');

%获取色度大和光度较弱的图像
[m,n,~] = size(I_Board);
I_UnW=double(zeros(m,n,3));
for i=1:m
    for j=1:n
        if((I_Board_hsv(i,j,2)>0.45 && I_Board_hsv(i,j,1)<0.6)==1)
            I_UnW(i,j,1)=I_Board_hsv(i,j,1);
            I_UnW(i,j,2)=I_Board_hsv(i,j,2);
            I_UnW(i,j,3)=I_Board_hsv(i,j,3)+0.4;
        end
    end
end

I_UnW_rgb=hsv2rgb(I_UnW);
I_UnW_rgb=uint8(255.*I_UnW_rgb);

subplot(222),imshow(I_UnW_rgb),title('非白色胶囊板图像');

I_UnW_Gray=im2gray(I_UnW_rgb);
I_UnW_Gray=medfilt2(I_UnW_Gray,[3,3]);

se=ones(5,5);
I_UnW_Gray_Open=imopen(I_UnW_Gray,se);

%使用I_UnW_BW_maker作为maker，对I_UnW_BW_erode进行重构图像
I_UnW_Gray=imreconstruct(I_UnW_Gray_Open,I_UnW_Gray);
subplot(223),imshow(I_UnW_Gray),title('灰度图片');

I_UnW_BW=imbinarize(I_UnW_Gray,0.001);
subplot(224),imshow(I_UnW_BW),title('二值化后的灰度图片');

se=ones(8,8);
I_UnW_Gray_Close=imclose(I_UnW_BW,se);
figure(5),subplot(221),imshow(I_UnW_Gray_Close),title('闭运算后的灰度图片');

se=ones(16,16);
I_UnW_Gray_Erode=imerode(I_UnW_Gray_Close,se);
subplot(222),imshow(I_UnW_Gray_Erode),title('腐蚀后的灰度图片');

%使用I_UnW_BW_maker作为maker，对I_UnW_BW_erode进行重构图像
I_UnW_Gray=imreconstruct(I_UnW_Gray_Erode,I_UnW_Gray_Close);
subplot(223),imshow(I_UnW_Gray),title('腐蚀后的灰度图片');

%统计标注连通域
%使用外接矩形框选连通域，并使用形心确定连通域位置
[l,b_n] = bwlabel(I_UnW_Gray,8);

status_box=regionprops(l,'BoundingBox');
centroid_box = regionprops(l,'Centroid');

%%重新排序
for i=1:b_n
    centriod(i,1:2)=centroid_box(i,1).Centroid(1,:);
    centriod(i,3)=i;
%     centriod(i,4)=0;
end

b_a=0;
for b=1:b_n
    centriod_min=min(centriod(:,2));
    for i=1:b_n
        if centriod_min > centriod(i,2)
            centriod_min=centriod(i,2);
        end
    end

    flag=1;
    for i=1:b_n
        if (centriod(i,2) > (centriod_min - 15)) && (centriod(i,2) < (centriod_min + 50))
            centriod_r(flag,:,b)=centriod(i,:);
            centriod(i,1)=1000;
            centriod(i,2)=1000;
            b_a=b_a+1;
            flag=flag+1;
        end
    end
    
    for i=1:flag-1
        for j=1:flag-1-i
            if(centriod_r(j,1,b)>centriod_r(j+1,1,b))
                [centriod_r(j,:,b),centriod_r(j+1,:,b)]=swap(centriod_r(j,:,b),centriod_r(j+1,:,b));
            end
        end
    end
    
    if b_a==b_n
        break;
    end
end

% template=load('set2template3.mat');
template=load('template.mat');
template=template.centriod_r;

[m,n,p]=size(centriod_r);
centriod_r_r=zeros(size(template));
q=1;
for i=1:p
    flag=0;
    for j=1:m
        if(centriod_r(j,1,i))
            flag=flag+1;
        end
    end
    if flag>=2
        for n=1:m
            centriod_r_r(n,:,q)=centriod_r(n,:,i);
        end
        q=q+1;
    end
end

[m,n,p]=size(template);
centriod_re=template;
for i=1:p
    for j=1:m
        centriod_re(j,3,i)=100;
    end
end
for i=1:p
    for j=1:m
        for o=1:m
            if (centriod_r_r(j,1,i) > (template(o,1,i) - 40)) && (centriod_r_r(j,1,i) < (template(o,1,i) + 30))
                centriod_re(o,:,i)=centriod_r_r(j,:,i);
            end 
        end
    end
end

[m,n,p]=size(centriod_r_r);
centriod_t=centriod;
for i=1:p
    for j=1:m
        centriod((i-1)*m+j,:)=centriod_re(j,:,i);
        centriod_t((i-1)*m+j,:)=template(j,:,i);
    end
end

[b_t,~]=size(centriod);

template_box=load('template_box.mat');
% template_box=load('set2template_box');
template_box=template_box.status_box;

subplot(224),imshow(I_Board_nor),title('板标记图片');hold on;

for i=1:b_t
    if centriod(i,3)~=100
        if status_box(centriod(i,3)).BoundingBox(1,3)>template_box(centriod(i,3)).BoundingBox(1,3)-25 && status_box(centriod(i,3)).BoundingBox(1,3)<template_box(centriod(i,3)).BoundingBox(1,3)+40 && status_box(centriod(i,3)).BoundingBox(1,4)>template_box(centriod(i,3)).BoundingBox(1,4)-17 && status_box(centriod(i,3)).BoundingBox(1,4)<template_box(centriod(i,3)).BoundingBox(1,4)+40
            rectangle('position',status_box(centriod(i,3)).BoundingBox,'edgecolor','g');
        else
            rectangle('position',template_box(centriod_t(i,3)).BoundingBox,'edgecolor','b');
%             centriod_t(i,3)
            template_box(centriod_t(i,3)).BoundingBox
        end
    else
        rectangle('position',template_box(centriod_t(i,3)).BoundingBox,'edgecolor','r');
    end
end

BoundingBox=zeros(1,2);
for i=1:b_n
    if status_box(i).BoundingBox(1,3)>BoundingBox(1,1)
        BoundingBox(1,1)=status_box(i).BoundingBox(1,3)+5;
    end
    if status_box(i).BoundingBox(1,4)>BoundingBox(1,2)
        BoundingBox(1,2)=status_box(i).BoundingBox(1,4)+5;
    end
end
BoundingBox=uint16(BoundingBox);

se=ones(5,11);
I_UnW_Gray_Close=imopen(I_UnW_Gray_Close,se);
for i=1:b_t
    if centriod(i,3)~=100
        if status_box(centriod(i,3)).BoundingBox(1,3)>template_box(centriod(i,3)).BoundingBox(1,3)-25 && status_box(centriod(i,3)).BoundingBox(1,3)<template_box(centriod(i,3)).BoundingBox(1,3)+15 && status_box(centriod(i,3)).BoundingBox(1,4)>template_box(centriod(i,3)).BoundingBox(1,4)-17 && status_box(centriod(i,3)).BoundingBox(1,4)<template_box(centriod(i,3)).BoundingBox(1,4)+15
            temp=I_UnW_Gray_Close(uint16(status_box(centriod(i,3)).BoundingBox(1,2)):uint16(status_box(centriod(i,3)).BoundingBox(1,2)+status_box(centriod(i,3)).BoundingBox(1,4)),uint16(status_box(centriod(i,3)).BoundingBox(1,1)):uint16(status_box(centriod(i,3)).BoundingBox(1,1)+status_box(centriod(i,3)).BoundingBox(1,3)),:);
        else
            temp=I_UnW_Gray_Close(uint16(template_box(centriod_t(i,3)).BoundingBox(1,2)):uint16(template_box(centriod_t(i,3)).BoundingBox(1,2)+template_box(centriod_t(i,3)).BoundingBox(1,4)),uint16(template_box(centriod_t(i,3)).BoundingBox(1,1)):uint16(template_box(centriod_t(i,3)).BoundingBox(1,1)+template_box(centriod_t(i,3)).BoundingBox(1,3)),:);
        end
    else
        temp=I_UnW_Gray_Close(uint16(template_box(centriod_t(i,3)).BoundingBox(1,2)):uint16(template_box(centriod_t(i,3)).BoundingBox(1,2)+template_box(centriod_t(i,3)).BoundingBox(1,4)),uint16(template_box(centriod_t(i,3)).BoundingBox(1,1)):uint16(template_box(centriod_t(i,3)).BoundingBox(1,1)+template_box(centriod_t(i,3)).BoundingBox(1,3)),:);
    end
%     temp=imresize(temp,BoundingBox);
    I_cell(i).cell(:,:,:)=temp;
end

%%连通域面积检测
figure(2);
[m,n,p]=size(centriod_r_r);
area=zeros(b_t,1);
for i=1:b_t
    subplot(p,m,i),imshow(I_cell(i).cell);
    Area_temp=0;
    [mm,nn]=size(I_cell(i).cell);
    for ii=1:mm
       for jj=1:nn
           if I_cell(i).cell(ii,jj)
               Area_temp=Area_temp+1;
           end
       end
    end

    area(i,1)=Area_temp;
    [r,c]=size(I_cell(i).cell);
    area(i,1)=area(i,1)/(r*c);
end

defect_flag=zeros(p,m);
for t=1:p
    for i=1:m
        if(area(i+(t-1)*m,1)<0.39)
           defect_flag(t,i)=1;
        end
    end
end

%%RGB直方图检测
%获取原始图片
for i=1:b_t
    if centriod(i,3)~=100
        if status_box(centriod(i,3)).BoundingBox(1,3)>template_box(centriod(i,3)).BoundingBox(1,3)-25 && status_box(centriod(i,3)).BoundingBox(1,3)<template_box(centriod(i,3)).BoundingBox(1,3)+15 && status_box(centriod(i,3)).BoundingBox(1,4)>template_box(centriod(i,3)).BoundingBox(1,4)-17 && status_box(centriod(i,3)).BoundingBox(1,4)<template_box(centriod(i,3)).BoundingBox(1,4)+15
            temp=I_Board_nor(uint16(status_box(centriod(i,3)).BoundingBox(1,2)):uint16(status_box(centriod(i,3)).BoundingBox(1,2)+status_box(centriod(i,3)).BoundingBox(1,4)),uint16(status_box(centriod(i,3)).BoundingBox(1,1)):uint16(status_box(centriod(i,3)).BoundingBox(1,1)+status_box(centriod(i,3)).BoundingBox(1,3)),:);
        else
            temp=I_Board_nor(uint16(template_box(centriod_t(i,3)).BoundingBox(1,2)):uint16(template_box(centriod_t(i,3)).BoundingBox(1,2)+template_box(centriod_t(i,3)).BoundingBox(1,4)),uint16(template_box(centriod_t(i,3)).BoundingBox(1,1)):uint16(template_box(centriod_t(i,3)).BoundingBox(1,1)+template_box(centriod_t(i,3)).BoundingBox(1,3)),:);
        end
    else
        temp=I_Board_nor(uint16(template_box(centriod_t(i,3)).BoundingBox(1,2)):uint16(template_box(centriod_t(i,3)).BoundingBox(1,2)+template_box(centriod_t(i,3)).BoundingBox(1,4)),uint16(template_box(centriod_t(i,3)).BoundingBox(1,1)):uint16(template_box(centriod_t(i,3)).BoundingBox(1,1)+template_box(centriod_t(i,3)).BoundingBox(1,3)),:);
    end
%     temp=imresize(temp,BoundingBox);
    Src_cell(i).cell(:,:,:)=temp;
end

figure(3);
[m,n,p]=size(centriod_r_r);
for i=1:b_t
    subplot(p,m,i),imshow(Src_cell(i).cell);
end
%%RGB 直方图检测
for t=1:p
   for i=1:m
       for j=i:m
           similarity_para(i,j,t)=same_detect(Src_cell(i+(t-1)*m).cell,Src_cell(j+(t-1)*m).cell,1);
       end
   end
   similarity_para(:,:,t)=similarity_para(:,:,t)-eye(m,m);
   similarity_para(:,:,t)=similarity_para(:,:,t)+similarity_para(:,:,t)';
end

[m,n,p]=size(similarity_para);
for t=1:p
    sim_max=0;
    for i=1:m
       for j=1:n
           if sim_max<similarity_para(i,j,t)
              sim_max=similarity_para(i,j,t);
           end
       end
    end
    sim_max_para(t,1)=sim_max;
end

similarity_flag=zeros(m,n,p);
for t=1:p
    for i=1:m
       for j=1:n
           if ((sim_max_para(t,1)-0.2>similarity_para(i,j,t)) && (similarity_para(i,j,t)~=0) || (similarity_para(i,j,t)<0.4)) && ~(similarity_para(i,j,t)>0.45)
              similarity_flag(i,j,t)=1;
           end
       end
    end
end

defect_rgb_flag=zeros(p,m);
for t=1:p
    for i=1:m
        temp_s=sum(similarity_flag(:,:,t),1);
        k=find(temp_s(:)>=(m));
        [r,~]=size(k);
        for j=1:r
            defect_rgb_flag(t,k(j))=1;
        end
    end
end

defect_mux_flag=defect_flag+defect_rgb_flag;

for t=1:p
    for i=1:m
        if(defect_mux_flag(t,i)>=1)
            defect_mux_flag(t,i)=1;
        else
            defect_mux_flag(t,i)=0;
        end
    end
end
%%矩阵重新排列
for t=1:p
    for i=1:m
        defect_flag_end((t-1)*m+i,1)=defect_mux_flag(t,i);
    end
end

figure(4),imshow(I),hold on;
status_box_b(1).BoundingBox(1,1)=10+status_box_b(1).BoundingBox(1,1);
status_box_b(1).BoundingBox(1,2)=25+status_box_b(1).BoundingBox(1,2);
status_box_b(1).BoundingBox(1,3)=0;
status_box_b(1).BoundingBox(1,4)=0;

for i=1:b_t
    if centriod(i,3)~=100
        if status_box(centriod(i,3)).BoundingBox(1,3)>template_box(centriod(i,3)).BoundingBox(1,3)-25 && status_box(centriod(i,3)).BoundingBox(1,3)<template_box(centriod(i,3)).BoundingBox(1,3)+40 && status_box(centriod(i,3)).BoundingBox(1,4)>template_box(centriod(i,3)).BoundingBox(1,4)-17 && status_box(centriod(i,3)).BoundingBox(1,4)<template_box(centriod(i,3)).BoundingBox(1,4)+40
            printlabel(status_box(centriod(i,3)).BoundingBox + status_box_b(1).BoundingBox,defect_flag_end(i));
        else
            printlabel(template_box(centriod_t(i,3)).BoundingBox + status_box_b(1).BoundingBox,defect_flag_end(i));
        end
    else
        printlabel(template_box(centriod_t(i,3)).BoundingBox + status_box_b(1).BoundingBox,defect_flag_end(i));
    end
end
clear all;
end