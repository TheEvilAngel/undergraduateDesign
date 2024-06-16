function FullFileNameArray=read_Feko_ffe_data(num_sam,flag) 
%#########################################################################
% Description: read  the *.ffe  files of Feko 

%Feko 默认的数据输出格式为.out文件，里面包含大量与计算数据无关的信息
%在FarField中可以设置，将Feko的数据输出到ffe文件中

%FEKO .ffe文件中远场电磁散射数据的排列方式
%**********文件**********
%**********频率**********
%**********极化**********
%*****电磁波入射角度*****
%*****电磁波散射角度*****

% Input parameters: 
%path-数据保存的路径
% num_sam - the number samples in elevation,azimuth,frequency,and files:  [Ma*Ne,Pa*Qe,Pol,Nf,NFile]
%flag: flag=1单基地散射 flag=2双基地散射
%Ma:入射电磁波扫描的方位角个数
%Ne:入射电磁波扫描的俯仰角个数
%Pa:散射电磁波扫描的方位角个数
%Qe:散射电磁波扫描的俯仰角个数
%Pol:电磁波入射极化通道，Pol=1:VH,HH或HV,VV；Pol=2:全极化，HH,VH,HV,VV
%Nf:频率扫描的点数
%NFile:ffe文件的个数

%单基地情况，收发一致，则Ma*Ne=Pa*Qe
%双基地情况，收发分置，则Ma*Ne~=Pa*Qe或Ma*Ne=Pa*Qe

%##########################################################################

[File,PathName]=uigetfile('*.ffe', 'Select the FEKO data file to open','MultiSelect','on'); %对话框选择多个文件
% PathName=uigetdir; %对话框选择指定文件
% PathName=strcat(PathName,'\');
Ne = num_sam(1);          % 电磁波入射方向角度个数
Na = num_sam(2);          % 电磁波散射方向角度个数
Pol=num_sam(3);           %入射电磁波极化通道个数
Nf = num_sam(4);          % 电磁波频率点数
NFile= num_sam(5);        % ffe.文件个数——适用于Matlab与Feko混合编程情况

if NFile~=1
    if length(File)~=NFile
        error('选中的文件个数与设置不一致')
    end
end

if flag==1
    N=Ne;   %单基地情况：电磁波入射方向和散射方向一致，电磁波角度维度为2维
    MatDescribe='后向散射-远场电磁计算数据';
    DirectIncident=strcat('电磁波入射角度个数:',num2str(Ne));
    DirectReceive=strcat('电磁波散射角度个数:',num2str(Na));
    Polarization=strcat('入射极化通道个数:',num2str(Pol));
    FrequencyPoint=strcat('频率扫描点数:',num2str(Nf));
elseif flag==2
    N=Ne*Na;%双基地情况：电磁波入射方向和散射方向不一致，电磁波角度维度为4维
    MatDescribe='双基地散射-电磁计算数据';
    DirectIncident=strcat('电磁波入射角度个数:',num2str(Ne));
    DirectReceive=strcat('电磁波散射角度个数:',num2str(Na));
    Polarization=strcat('入射极化通道个数:',num2str(Pol));
    FrequencyPoint=strcat('频率扫描点数:',num2str(Nf));
end

if Pol==1
        data=cell(2+5,1);
elseif Pol==2
        data=cell(4+5,1);
end
    
h=waitbar(0,'please wait');

for ProjectNum=1:NFile
     if NFile~=1
        FileName=File{ProjectNum};
    else
        FileName=File;
     end
    
    fid = fopen([PathName,FileName], 'r'); 
    
    for j=1:Nf
        for p=1:Pol
            for i=1:N
                tLine=fgetl(fid);  % 读取第一行数据
                while contains(tLine,'#')||contains(tLine,'*')||isempty(tLine)%数据行是不包含符号'#'、'*'，并且也不可能是空行
                    biteNull=ftell(fid);%得到指针当前位置相对与文件首的偏移字节数
                    tLine=fgetl(fid);  % 继续读取下一行数据
                end
                fseek(fid,biteNull,'bof');%定位到第一行数据的位置
                %数据为9列：Theta,Phi,Re(Etheta),Im(Etheta),Re(Ephi),Im(Ephi),RCS(Theta),RCS(Phi),RCS(Total)
                A = fscanf(fid,'%f %f %f %f %f %f %f %f %f',[1 9]);

                E_real_H = A(5);% H
                E_imag_H = A(6);
                E_real_V = A(3);% V
                E_imag_V = A(4);

                RCS_complex_H(i,p,j) = E_real_H+1j*E_imag_H;
                RCS_complex_V(i,p,j) = E_real_V+1j*E_imag_V;
                
                display=['运行中...',num2str((Nf*Pol*N*(ProjectNum-1)+(Pol*N*(j-1)+(N*(p-1)+i)))/NFile/Nf/Pol/N*100),'%'];%进度条显示
                waitbar((Nf*Pol*N*(ProjectNum-1)+(Pol*N*(j-1)+(N*(p-1)+i)))/NFile/Nf/Pol/N,h,display)
                
            end
        end
    end
    
    %数据分割及保存
    if Pol==1
        if flag==1
            data{5+1,1}=squeeze(RCS_complex_H).';
            data{5+2,1}=squeeze(RCS_complex_V).';
        else
            data{5+1,1}=permute(reshape(squeeze(RCS_complex_H),Na,Ne,[]),[2,1,3]);
            data{5+2,1}=permute(reshape(squeeze(RCS_complex_V),Na,Ne,[]),[2,1,3]);
        end
    elseif Pol==2
        data{5+1,1}=permute(reshape(squeeze(RCS_complex_H(:,1,:)),Na,Ne,[]),[2,1,3]);
        data{5+2,1}=permute(reshape(squeeze(RCS_complex_V(:,1,:)),Na,Ne,[]),[2,1,3]);
        data{5+3,1}=permute(reshape(squeeze(RCS_complex_H(:,2,:)),Na,Ne,[]),[2,1,3]);
        data{5+4,1}=permute(reshape(squeeze(RCS_complex_V(:,2,:)),Na,Ne,[]),[2,1,3]);
    end
    
    data{1,1}=MatDescribe;
    data{2,1}=DirectIncident;
    data{3,1}=DirectReceive;
    data{4,1}=Polarization;
    data{5,1}=FrequencyPoint;
    
    FileNameFrontCell=strsplit(FileName,'.');
    FileNameFront= FileNameFrontCell{1};%文件名前缀
    FullFileName=strcat(PathName,FileNameFront,'.mat');
    save(FullFileName,'data');%按文件将数据进行保存
    FullFileNameArray{ProjectNum}= FullFileName;
    clearvars RCS_complex_H  RCS_complex_V data%删除数据，防止内存不够
    fclose(fid); %每个文件读完一定要记得马上关闭
end
end

 


