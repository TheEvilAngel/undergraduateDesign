function FullFileNameArray=read_Feko_ffe_data(num_sam,flag) 
%#########################################################################
% Description: read  the *.ffe  files of Feko 

%Feko Ĭ�ϵ����������ʽΪ.out�ļ������������������������޹ص���Ϣ
%��FarField�п������ã���Feko�����������ffe�ļ���

%FEKO .ffe�ļ���Զ�����ɢ�����ݵ����з�ʽ
%**********�ļ�**********
%**********Ƶ��**********
%**********����**********
%*****��Ų�����Ƕ�*****
%*****��Ų�ɢ��Ƕ�*****

% Input parameters: 
%path-���ݱ����·��
% num_sam - the number samples in elevation,azimuth,frequency,and files:  [Ma*Ne,Pa*Qe,Pol,Nf,NFile]
%flag: flag=1������ɢ�� flag=2˫����ɢ��
%Ma:�����Ų�ɨ��ķ�λ�Ǹ���
%Ne:�����Ų�ɨ��ĸ����Ǹ���
%Pa:ɢ���Ų�ɨ��ķ�λ�Ǹ���
%Qe:ɢ���Ų�ɨ��ĸ����Ǹ���
%Pol:��Ų����伫��ͨ����Pol=1:VH,HH��HV,VV��Pol=2:ȫ������HH,VH,HV,VV
%Nf:Ƶ��ɨ��ĵ���
%NFile:ffe�ļ��ĸ���

%������������շ�һ�£���Ma*Ne=Pa*Qe
%˫����������շ����ã���Ma*Ne~=Pa*Qe��Ma*Ne=Pa*Qe

%##########################################################################

[File,PathName]=uigetfile('*.ffe', 'Select the FEKO data file to open','MultiSelect','on'); %�Ի���ѡ�����ļ�
% PathName=uigetdir; %�Ի���ѡ��ָ���ļ�
% PathName=strcat(PathName,'\');
Ne = num_sam(1);          % ��Ų����䷽��Ƕȸ���
Na = num_sam(2);          % ��Ų�ɢ�䷽��Ƕȸ���
Pol=num_sam(3);           %�����Ų�����ͨ������
Nf = num_sam(4);          % ��Ų�Ƶ�ʵ���
NFile= num_sam(5);        % ffe.�ļ���������������Matlab��Feko��ϱ�����

if NFile~=1
    if length(File)~=NFile
        error('ѡ�е��ļ����������ò�һ��')
    end
end

if flag==1
    N=Ne;   %�������������Ų����䷽���ɢ�䷽��һ�£���Ų��Ƕ�ά��Ϊ2ά
    MatDescribe='����ɢ��-Զ����ż�������';
    DirectIncident=strcat('��Ų�����Ƕȸ���:',num2str(Ne));
    DirectReceive=strcat('��Ų�ɢ��Ƕȸ���:',num2str(Na));
    Polarization=strcat('���伫��ͨ������:',num2str(Pol));
    FrequencyPoint=strcat('Ƶ��ɨ�����:',num2str(Nf));
elseif flag==2
    N=Ne*Na;%˫�����������Ų����䷽���ɢ�䷽��һ�£���Ų��Ƕ�ά��Ϊ4ά
    MatDescribe='˫����ɢ��-��ż�������';
    DirectIncident=strcat('��Ų�����Ƕȸ���:',num2str(Ne));
    DirectReceive=strcat('��Ų�ɢ��Ƕȸ���:',num2str(Na));
    Polarization=strcat('���伫��ͨ������:',num2str(Pol));
    FrequencyPoint=strcat('Ƶ��ɨ�����:',num2str(Nf));
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
                tLine=fgetl(fid);  % ��ȡ��һ������
                while contains(tLine,'#')||contains(tLine,'*')||isempty(tLine)%�������ǲ���������'#'��'*'������Ҳ�������ǿ���
                    biteNull=ftell(fid);%�õ�ָ�뵱ǰλ��������ļ��׵�ƫ���ֽ���
                    tLine=fgetl(fid);  % ������ȡ��һ������
                end
                fseek(fid,biteNull,'bof');%��λ����һ�����ݵ�λ��
                %����Ϊ9�У�Theta,Phi,Re(Etheta),Im(Etheta),Re(Ephi),Im(Ephi),RCS(Theta),RCS(Phi),RCS(Total)
                A = fscanf(fid,'%f %f %f %f %f %f %f %f %f',[1 9]);

                E_real_H = A(5);% H
                E_imag_H = A(6);
                E_real_V = A(3);% V
                E_imag_V = A(4);

                RCS_complex_H(i,p,j) = E_real_H+1j*E_imag_H;
                RCS_complex_V(i,p,j) = E_real_V+1j*E_imag_V;
                
                display=['������...',num2str((Nf*Pol*N*(ProjectNum-1)+(Pol*N*(j-1)+(N*(p-1)+i)))/NFile/Nf/Pol/N*100),'%'];%��������ʾ
                waitbar((Nf*Pol*N*(ProjectNum-1)+(Pol*N*(j-1)+(N*(p-1)+i)))/NFile/Nf/Pol/N,h,display)
                
            end
        end
    end
    
    %���ݷָ����
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
    FileNameFront= FileNameFrontCell{1};%�ļ���ǰ׺
    FullFileName=strcat(PathName,FileNameFront,'.mat');
    save(FullFileName,'data');%���ļ������ݽ��б���
    FullFileNameArray{ProjectNum}= FullFileName;
    clearvars RCS_complex_H  RCS_complex_V data%ɾ�����ݣ���ֹ�ڴ治��
    fclose(fid); %ÿ���ļ�����һ��Ҫ�ǵ����Ϲر�
end
end

 


