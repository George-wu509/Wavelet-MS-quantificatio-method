 function Peak_IDQ_3(mat,x,C)
%ex: Peak_IDQ_3(0,0,0);        for inner simulated peaks data (choose1) 
%ex: Peak_IDQ_3(1,0,0);        To load all .mat file in dictionary /data source (choose2) 
%ex: Peak_IDQ_3(xxx.mat,0,0);  To load xxx.mat in dictionary /data source (choose3) 
%ex: Peak_IDQ_3(0,x,0);        for existed variable x in workspace (choose4) 
%ex: Peak_IDQ_3(0,0,C);        for existed variable C in workspace (choose5) 
% Edit: 2012 0802 am11  by Ching

if exist(strcat(cd,'Coding sub'),'dir')==0
    addpath('Coding sub');
end

%% 0 ===== parameter =====

% 0_1 -- make peak parameters -- 
func=1;                   %func=1 gaussian function  func=2  cachy function  func=3 mixed gaussian and cachy functions
dataLength=5000;
peakN=15;
w_bound=30;
h_bound=1000;
noisy_value=15;
repeat_time=3;

% 0_2 -- cwt parameters -- 
startS=1;
endS=30;
s_range_C=startS:endS;
s_range_Area=5:30;

% 0_3 -- RidgeID parameters -- 
mm_window=0;              % mm_window: window size of min/max sliding window
SNR_window=1000;          % SNR_window: window size of SNR sliding window
N_adover=2;               % max error adjacent tolerant in over-scales analysis 0,1,2

% 0_4 -- batch parameters -- 
gap_batch=0.5; %0.3:0.2:0.9;            % gap: threhold of min/max value of CWT coefficient
gap_snr_batch=5; %3:2:7;          % gap_snr: cwt(s,b)/noise_level>gap_snr as meaningful peak

% 0_5 -- cwt Area parameters -- 
Area_adover=N_adover+1;     % Area_adover: Calculate cwt Area using parameter adover
overlap=0;                % overlap=0: Neglect overlap effect, overlap=1: calculate with overlap effect
ratio_matrix_n=1;         % ratio_matrix_n=1: load Gauss.mat,  ratio_matrix_n=2: load Cauchy.mat, ratio_matrix_n=3: load GauCauchy.mat
IDsum_gap=5;             % IDsum_gap: threadhole in IDsum which represented Peaks 
hat_ratio=5/1.7;          % wavelet function ratio of "center to min" : "center to end of function" 
C_ov_tempL=500;           % C_ov_tempL: temp of C matrix when calculate overlap peaks.
pmin_number=10;

%% 1 ===== Main Program =====
if mat==1
    choose=2;
elseif mat==0&&isscalar(x)&&isscalar(C)
    choose=1;
else
    show('error');return;
end
if ratio_matrix_n==1
    eval(['load(''','data result/','Gauss.mat',''');']);
elseif ratio_matrix_n==2
    eval(['load(''','data result/','Cauchy.mat',''');']);
elseif ratio_matrix_n==3
    eval(['load(''','data result/','GauCauchy.mat',''');']);
end
        
%% 11 ===== Peak_IDQ_1(0,0,0)  for inner simulated peaks data =====
if choose==1 %ex: Peak_IDQ_1(0,0,0);        for inner simulated peaks data (choose1)    
    t_c=0;
    for re=1:repeat_time        
        % 11_1 -- IO  --
        clear h1 h2 h3 h4 ID Ridge_ID Valley_ID IDC IDN Ridge x C;
        cl=clock;
        new_fdir0=strcat('data result/');
        new_fdir1=strcat(int2str(cl(2)),int2str(cl(3)),'_',int2str(cl(4)),'_',int2str(cl(5)),'_no',int2str(re),'/');
        new_fdir2=[new_fdir0,new_fdir1];
        
        % 11_2 -- make peak  --
        xx=0;ParaM=0;DrawPeak=0;
        if dataLength<=SNR_window
            SNR_window=dataLength;
        end
        [x,ParaM,DrawPeak]=Randmake_peak(dataLength,peakN,w_bound,h_bound,xx,ParaM,DrawPeak,func);                    
        xn=[x(1,:);awgn(x(2,:),noisy_value,'measured');x(3,:)];

        % 11_3 -- cwt -- 
        C=cwt(x(2,1:dataLength),s_range_C,'mexh');
        Cn=cwt(xn(2,1:dataLength),s_range_C,'mexh');
        for pn=1:peakN; 
            Cisop{pn}=cwt(DrawPeak(pn,1:dataLength),s_range_C,'mexh');
        end
        t_count=length(gap_batch)*length(gap_snr_batch)*repeat_time;
        
        % 11_4 -- bacth run and ridge attract  --
        for gap=gap_batch  %0.3:0.2:0.9;            % gap: threhold of min/max value of CWT coefficient
            for gap_snr=gap_snr_batch  %3:2:7;          % gap_snr: cwt(s,b)/noise_level>gap_snr as meaningful peak
                tic;
                t_c=t_c+1;
                % 11_4_1 -- IO  --
                new_fdir3=strcat('gap_',int2str(gap*10),'_snrG_',int2str(gap_snr*10),'/');
                new_fdir=[new_fdir2,new_fdir3];
                if exist(new_fdir,'dir')==0
                    mkdir(new_fdir);
                end
                new_fname=strcat('N',num2str(peakN),'L',num2str(dataLength),'.mat');                
                disp(strcat('Process ( ',int2str(t_c),'/',int2str(t_count),')'));
                disp(strcat('gap:  ',mat2str(gap*10),'    gap_snr:  ',mat2str(gap_snr*10)));
                
                % 11_4_2 -- ridge attract  --
                [ID,IDW,IDsum,snr,noisy_level]=RidgeID(C,startS,endS,mm_window,SNR_window,gap_snr,gap,N_adover);
                [ratio_matrix,ratio,cwt_Area,ParaM_p,IDW_ov,IDsum_ov,C_ov]=CWT_Area(ID,IDW{Area_adover},IDsum{Area_adover},C,D_dist,C50,Area,s_range,s_range_Area,w_range,IDsum_gap,hat_ratio,C_ov_tempL,SNR_window,gap_snr,gap,N_adover,overlap,ratio_matrix_n,pmin_number);
                [IDn,IDWn,IDsumn,snrn,noisy_leveln]=RidgeID(Cn,startS,endS,mm_window,SNR_window,gap_snr,gap,N_adover);
                [ratio_matrixn,ration,cwt_Arean,ParaM_pn,IDW_ovn,IDsum_ovn,C_ovn]=CWT_Area(IDn,IDWn{Area_adover},IDsumn{Area_adover},Cn,D_dist,C50,Area,s_range,s_range_Area,w_range,IDsum_gap,hat_ratio,C_ov_tempL,SNR_window,gap_snr,gap,N_adover,overlap,ratio_matrix_n,pmin_number);
                
                % 11_4_3 -- Save files  --
                eval(['save(''',new_fdir,new_fname,''',''x'',''xn'',''ParaM'',''DrawPeak'',''C'',''Cn'',''Cisop'');']);
                eval(['save(''',new_fdir,'IDresult',''',''ID'',''IDn'',''IDW'',''IDWn'',''IDsum'',''IDsumn'',''snr'',''snrn'',''noisy_level'',''noisy_leveln'');']);
                eval(['save(''',new_fdir,'Area',''',''ratio_matrix'',''ratio_matrixn'',''ratio'',''ration'',''cwt_Area'',''cwt_Arean'',''ParaM_p'',''ParaM_pn'');']);
                eval(['save(''',new_fdir,'Check_overlap',''',''C_ov'',''IDsum_ov'',''IDW_ov'');']);
                toc;
            end
        end
    end %end of re
%% 12 ===== Peak_IDQ_1(1,0,0) Dictionary batch run in data source =====
elseif choose==2  %Peak_IDQ_1(1,0,0);        To load all .mat file in dictionary /data source (choose2) 
    
    % 12_1 -- IO and list  --
    new_fdir0=strcat('data result/');
    new_fdir2=new_fdir0;
    source0=dir('data source/*.mat');
    [a,a0]=size(source0);
    source_files={source0(:,1).name};clear source0 a0;
    
    % 12_2 -- list run for every *.mat  --
    t_count=(length(gap_batch))*(length(gap_snr_batch))*a;t_c=0;
    for diri=1:a    % each .mat files in /data source sictionary.
        disp(strcat('dataset:  ',source_files{1,diri}));               
                
        % 12_2_1 -- Read dictionary information  --
        clear h1 h2 h3 h4 ID Ridge_ID Valley_ID IDC IDN Ridge x C;
        eval(['matV=who(''-file'',''',cd,'/data source/',source_files{1,diri},''');']);
        eval(['load(''','data source/',source_files{1,diri},''');']);
        new_fname=strrep(source_files{1,diri},'.mat','_PI.mat');
        new_fname2=strrep(source_files{1,diri},'.mat','_IDresult');
        new_fname3=strrep(source_files{1,diri},'.mat','/');
        new_fname4=strrep(source_files{1,diri},'.mat','_Area');
        if exist('C','var')
        else
            C=cwt(x(2,:),s_range_C,'mexh');
        end
        for gap=gap_batch
            for gap_snr=gap_snr_batch
                t_c=t_c+1;tic;
                disp(strcat('Process ( ',int2str(t_c),'/',int2str(t_count),')'));
                new_fdir=[new_fdir2,'DIR/',new_fname3,strcat('gap_',int2str(gap*10),'_snrG_',int2str(gap_snr*10),'/')];
                if exist(new_fdir,'dir')==0
                    mkdir(new_fdir);
                end
                eval(['copyfile(''','data source/',source_files{1,diri},''',''',new_fdir,new_fname,''')']);
                disp(strcat('gap:  ',mat2str(gap),'    gap_snr:  ',mat2str(gap_snr)));           
                [ID,IDW,IDsum,snr,noisy_level]=RidgeID(C,startS,endS,mm_window,SNR_window,gap_snr,gap,N_adover);
                eval(['save(''',new_fdir,new_fname2,''',''ID'',''C'',''IDW'',''IDsum'',''snr'',''noisy_level'');']);
                [ratio_matrix,ratio,cwt_Area,ParaM_p,IDW_ov,IDsum_ov,C_ov]=CWT_Area(ID,IDW{Area_adover},IDsum{Area_adover},C,D_dist,C50,Area,s_range,s_range_Area,w_range,IDsum_gap,hat_ratio,C_ov_tempL,SNR_window,gap_snr,gap,N_adover,overlap,ratio_matrix_n,pmin_number);
                eval(['save(''',new_fdir,new_fname4,''',''ratio_matrix'',''ratio'',''cwt_Area'',''ParaM_p'');']);
                toc;
            end
        end
    end % end of source file.mat
  
end % End of Choose

function [ID,IDW,IDsum,snr,noisy_level]=RidgeID(input,init_scale,end_scale,mm_window,SNR_window,gap_snr,gap,N_adover)

%% Method 0 (Environment)
[ma,mb]=size(input);
ID(1:N_adover+4)={zeros(end_scale-init_scale+1,mb)};
IDW(1:N_adover+1)={zeros(end_scale-init_scale+1,mb)};
IDsum(1:N_adover+1)={zeros(1,mb)};

%%  Method 1  (max/min method)
if exist(strcat(cd,'Coding sub'),'dir')==0
     addpath('Coding sub');
end
if mm_window==0
for s=init_scale:end_scale
    for i=s+1:mb-s
        [ia,ib]=min(input(s,i-s:i+s));
        if ib~=1&&ib~=2*s+1
            ID{end-2}(s-init_scale+1,ib-s+i-1)=-1;
        end
        [ia,ib]=max(input(s,i-s:i+s));
        if ib~=1&&ib~=2*s+1
            ID{end-2}(s-init_scale+1,ib-s+i-1)=1;
        end
    end
end
else
    for s=init_scale:end_scale
        for i=1:mb
            if i<mm_window/2
                [ia,ib]=min(input(s,1:i+mm_window/2-1));
                if ib~=1&&ib~=i+mm_window/2-1
                    ID{end-2}(s-init_scale+1,ib)=-1;
                end
                [ia,ib]=max(input(s,1:i+mm_window/2-1));
                if ib~=1&&ib~=i+mm_window/2-1
                    ID{end-2}(s-init_scale+1,ib)=1;
                end
            elseif i>mb-mm_window/2+1
                [ia,ib]=min(input(s,i-mm_window/2:mb));
                if ib~=1&&ib~=mb-i+mm_window/2+1
                    ID{end-2}(s-init_scale+1,ib-mm_window/2+i-1)=-1;
                end
                [ia,ib]=max(input(s,i-mm_window/2:mb));
                if ib~=1&&ib~=mb-i+mm_window/2+1
                    ID{end-2}(s-init_scale+1,ib-mm_window/2+i-1)=1;
                end
            else
                [ia,ib]=min(input(s,i-mm_window/2:i+mm_window/2-1));
                if ib~=1&&ib~=mm_window
                    ID{end-2}(s-init_scale+1,ib-mm_window/2+i-1)=-1;
                end
                [ia,ib]=max(input(s,i-mm_window/2:i+mm_window/2-1));
                if ib~=1&&ib~=mm_window
                    ID{end-2}(s-init_scale+1,ib-mm_window/2+i-1)=1;
                end
            end
        end
    end
end
%%  Method 2  (max/min method and SNR positive(95%))
[snr,noisy_level]=SNR(input,init_scale:end_scale,0,SNR_window,gap);
for s=init_scale:end_scale
    for i=1:mb
        if snr(s,i)>=gap_snr&&ID{end-2}(s-init_scale+1,i)==1
            ID{end-1}(s-init_scale+1,i)=1;
            [ID{end-1}(s-init_scale+1,:),ID{end}(s-init_scale+1,i)]=Wdist(ID{end-2}(s-init_scale+1,:),i,mb,ID{end-1}(s-init_scale+1,:));
        end
    end
end
%%  Method 3  (max/min method and SNR positive(95%)and negative(5%) and over scale and adjacent slope fix ) 
for adover=0:N_adover           % adover: adjacent tolerant in over-scales analysis
    IDcheck(1:2)={1};imax=zeros(1,2);smax=zeros(1,2);iiw(1:2)={0};
    for i=1:mb
        if mod(i,2000)==0
            disp(strcat('m/z data: (',int2str(i),' / ',int2str(mb),')'));
        end               
        id5=0;id5max=0;mbb=zeros(1,mb);ii(1)={i};
        for s=init_scale:end_scale
            if s==init_scale
                if any(ID{end-1}(1,max(1,i-adover):min(i+adover,mb))==1)==0
                    id5=0;mbb=zeros(1,mb);
                else
                    fresult=find(ID{end-1}(1,max(1,i-adover):min(i+adover,mb))==1); % fresult: in adover window in specific scale and i, positions of ID{5}=1
                    fresult=fresult+max(1,i-adover)-1;
                    [AB,logic]=match(fresult,i,adover);
                    if logic==0
                        id5=0;mbb=zeros(1,mb);ii{1}=i;
                    else
                        ii{1}=AB;
                        id5=id5+1;
                        for bb=1:length(ii{1})
                           mbb(1,ii{1}(bb))=mbb(1,ii{1}(bb))+1;                          
                        end
                        if id5>id5max
                           id5max=id5;
                           smax(1)=s;[mbb_a,imax(1)]=max(mbb);clear mbb_a;
                        end
                    end
                end
            else % s~=init_scale
                if any(ID{end-1}(s-init_scale+1,max(1,ii{s-init_scale}(1)-adover):min(ii{s-init_scale}(end)+adover,mb))==1)==0
                    id5=0;mbb=zeros(1,mb);ii{s-init_scale+1}=i;
                else
                    fresult=find(ID{end-1}(s-init_scale+1,max(1,ii{s-init_scale}(1)-adover):min(ii{s-init_scale}(end)+adover,mb))==1);
                    fresult=fresult+max(1,ii{s-init_scale}(1)-adover)-1;
                    [AB,logic]=match(fresult,ii{s-init_scale},adover);
                    if logic==0
                        id5=0;mbb=zeros(1,mb);ii{s-init_scale+1}=i;
                    else
                        ii{s-init_scale+1}=AB;
                        id5=id5+1;
                        for bb=1:length(ii{s-init_scale+1})
                           mbb(1,ii{s-init_scale+1}(bb))=mbb(1,ii{s-init_scale+1}(bb))+1;
                        end
                        if id5>id5max
                           id5max=id5;
                           smax(1)=s;[mbb_a,imax(1)]=max(mbb);clear mbb_a;
                        end
                    end
                 end                   
            end

        end % end of s loop           
        if id5max>0
            iiw{1}=0;
            for s=smax(1)-id5max+1:smax(1)
                [a,b]=min(abs(ii{s}-imax(1)));
                iiw{1}=[iiw{1} ID{end}(s,ii{s}(b))];
            end
            iiw{1}=iiw{1}(2:end);      
            if imax(2)==0
                IDcheck{1}=ii{smax(1)-init_scale+1};IDcheck{2}=id5max;smax(2)=smax(1);imax(2)=imax(1);iiw{2}=iiw{1};
            else                      
                [AB,logic]=match(IDcheck{1},ii{smax(1)-init_scale+1},adover);                                           
                if logic==0
                    IDsum{adover+1}(1,imax(2))=IDcheck{2};
                    ID{adover+1}(smax(2)-IDcheck{2}+1:smax(2),imax(2))=1;
                    if isempty(find(iiw{2},1))==1
                        IDW{adover+1}(smax(2)-IDcheck{2}+1:smax(2),imax(2))=0;
                    else
                        IDW{adover+1}(smax(2)-IDcheck{2}+1:smax(2),imax(2))=iiw{2}(1:IDcheck{2});
                    end
                    IDcheck{1}=ii{smax(1)-init_scale+1};IDcheck{2}=id5max;iiw{2}=iiw{1};smax(2)=smax(1);imax(2)=imax(1);
                elseif logic==1&&id5max>=IDcheck{2}
                    IDcheck{1}=ii{smax(1)-init_scale+1};IDcheck{2}=id5max;iiw{2}=iiw{1};smax(2)=smax(1);imax(2)=imax(1);
                end
            end % end /if imax(2)==0/              
        end % end /if id5max>0/
    end % end /i/
    if imax(2)==0||smax(2)==0
    else
        IDsum{adover+1}(1,imax(2))=IDcheck{2};
        ID{adover+1}(smax(2)-IDcheck{2}+1:smax(2),imax(2))=1;
        IDW{adover+1}(smax(2)-IDcheck{2}+1:smax(2),imax(2))=iiw{2}(1:IDcheck{2});
    end
end

function [ratio_matrix,ratio,cwt_Area,ParaM_p,IDW_ov,IDsum_ov,C_ov]=CWT_Area(ID,IDW,IDsum,C,D_dist,C50,Area,s_range,s_range_Area,w_range,IDsum_gap,hat_ratio,C_ov_tempL,SNR_window,gap_snr,gap,N_adover,overlap,ratio_matrix_n,pmin_number)

%% 01 Parameter Set
C_s1=1;        % C_s1: CWT scale init value in 2D-CWT coefficient matrix

%% 02 Def variables
IDWsum=zeros(1,length(IDW(1,:)));
for IDWsum_i=1:length(IDW(1,:))
    for s=1:length(IDW(:,1))
        if IDW(s,IDWsum_i)>0
            IDWsum(1,IDWsum_i)=IDWsum(1,IDWsum_i)+1;
        end
    end
end    
IDsum_b=find(IDWsum>=IDsum_gap);    % IDsum_b: Peaks location number in IDsum
s_range_C=C_s1:length(IDW(:,1)');  

%% 10 Main program
if overlap==0&&isempty(find(IDW,1))~=1
    [cwt_Area,ratio,ratio_matrix,ParaM_p]=CalArea(C,IDWsum,IDW,IDsum_b,D_dist,C50,Area,s_range,s_range_Area,w_range,pmin_number,ratio_matrix_n);
    IDW_ov=0;IDsum_ov=0;C_ov=0;
elseif overlap==1&&isempty(find(IDW,1))~=1
    C_ov=C;IDW_ov=IDW;IDsum_ov=IDWsum;
    IDsum_L=zeros(length(s_range_C),length(IDsum_b));
    IDsum_R=zeros(length(s_range_C),length(IDsum_b));
    IDsum_0=zeros(2,length(IDsum_b));
    for IDsum_b_i=1:length(IDsum_b)
        for IDsum_s=1:length(s_range_C)
            if IDW(IDsum_s,IDsum_b(IDsum_b_i))>0
                IDsum_L(IDsum_s:IDWsum(IDsum_b(IDsum_b_i)),IDsum_b_i)=max(IDsum_b(IDsum_b_i)-fix(IDW(IDsum_s:IDWsum(IDsum_b(IDsum_b_i)),IDsum_b(IDsum_b_i))/2*hat_ratio)+1,1);
                IDsum_R(IDsum_s:IDWsum(IDsum_b(IDsum_b_i)),IDsum_b_i)=min(IDsum_b(IDsum_b_i)+fix(IDW(IDsum_s:IDWsum(IDsum_b(IDsum_b_i)),IDsum_b(IDsum_b_i))/2*hat_ratio)+1,length(IDW));
                IDsum_0(1,IDsum_b_i)=IDsum_s;IDsum_0(2,IDsum_b_i)=IDsum(IDsum_b(IDsum_b_i));
                break;
            end
        end
    end

%% ------ [1] IDsum_b_over -----------------------------
 
    IDsum_b_over=zeros(2,1); % IDsum_b_over: An overlap peak set list (tempolar) 
    for IDsum_b_i=2:length(IDsum_b)
        if isempty(find(IDsum_L(max(IDsum_0(1,IDsum_b_i-1),IDsum_0(1,IDsum_b_i)):min(IDsum_0(2,IDsum_b_i-1),IDsum_0(2,IDsum_b_i)),IDsum_b_i)-IDsum_R(max(IDsum_0(1,IDsum_b_i-1),IDsum_0(1,IDsum_b_i)):min(IDsum_0(2,IDsum_b_i-1),IDsum_0(2,IDsum_b_i)),IDsum_b_i-1)<0,1))==0&&IDsum_b_i~=length(IDsum_b) % detect overlap
            IDsum_b_over=[IDsum_b_over [IDsum_b_i-1;IDsum_b(IDsum_b_i-1)]];
        else
            if length(IDsum_b_over(1,:))>1||(IDsum_b_i==length(IDsum_b)&&isempty(find(IDsum_L(max(IDsum_0(1,IDsum_b_i-1),IDsum_0(1,IDsum_b_i)):min(IDsum_0(2,IDsum_b_i-1),IDsum_0(2,IDsum_b_i)),IDsum_b_i)-IDsum_R(max(IDsum_0(1,IDsum_b_i-1),IDsum_0(1,IDsum_b_i)):min(IDsum_0(2,IDsum_b_i-1),IDsum_0(2,IDsum_b_i)),IDsum_b_i-1)<0,1))==0)
                if IDsum_b_i~=length(IDsum_b)||isempty(find(IDsum_L(max(IDsum_0(1,IDsum_b_i-1),IDsum_0(1,IDsum_b_i)):min(IDsum_0(2,IDsum_b_i-1),IDsum_0(2,IDsum_b_i)),IDsum_b_i)-IDsum_R(max(IDsum_0(1,IDsum_b_i-1),IDsum_0(1,IDsum_b_i)):min(IDsum_0(2,IDsum_b_i-1),IDsum_0(2,IDsum_b_i)),IDsum_b_i-1)<0,1))~=0
                    IDsum_b_over=[IDsum_b_over [IDsum_b_i-1;IDsum_b(IDsum_b_i-1)]];
                elseif isempty(find(IDsum_L(max(IDsum_0(1,IDsum_b_i-1),IDsum_0(1,IDsum_b_i)):min(IDsum_0(2,IDsum_b_i-1),IDsum_0(2,IDsum_b_i)),IDsum_b_i)-IDsum_R(max(IDsum_0(1,IDsum_b_i-1),IDsum_0(1,IDsum_b_i)):min(IDsum_0(2,IDsum_b_i-1),IDsum_0(2,IDsum_b_i)),IDsum_b_i-1)<0,1))==0
                    IDsum_b_over=[IDsum_b_over [IDsum_b_i-1;IDsum_b(IDsum_b_i-1)] [IDsum_b_i;IDsum_b(IDsum_b_i)]];                    
                end
                IDsum_b_over=IDsum_b_over(:,2:end);  % IDsum_b_over:  indenpedent overlap peaks in overlap datasets
                                                
%% --------[2] Analysis overlap peak sets in every scales-------
                if SNR_window>C_ov_tempL
                    SNR_window=C_ov_tempL;
                end                
                for i=1:fix(length(IDsum_b_over)/2)
                    if i==1
                        Ctemp=C;
                    else
                        Ctemp=Cother;
                    end
                    [Lpeak,Rpeak,w_L,w_R,Cother]=isolate(Ctemp,[IDsum_b_over(2,i) IDsum_b_over(2,end-i+1)],ID{end-1},ID{end},IDW,s_range_C,ratio_matrix_n,s_range_Area,pmin_number);
                    if ratio_matrix_n==1
                        eval(['load(''','data result/','Gauss.mat',''');']);
                    elseif ratio_matrix_n==2
                        eval(['load(''','data result/','Cauchy.mat',''');']);
                    elseif ratio_matrix_n==3
                        eval(['load(''','data result/','GauCauchy.mat',''');']);
                    end                   
                    C_ov(s_range_C,Lpeak(1,1))=C50(s_range_C,w_L)*Lpeak(3,1);C_ov(s_range_C,Rpeak(1,1))=C50(s_range_C,w_R)*Rpeak(3,1);
                    IDW_ov(s_range_C,Lpeak(1,1))=D_dist(s_range_C,w_L);IDW_ov(s_range_C,Rpeak(1,1))=D_dist(s_range_C,w_R);
                    IDsum_b(1,IDsum_b_over(1,i))=Lpeak(1,1);IDsum_b(1,IDsum_b_over(1,end-i+1))=Rpeak(1,1);
                end
                if rem(length(IDsum_b_over),2)~=0
                    C_ov(s_range_C,IDsum_b_over(2,((length(IDsum_b_over)+1)/2)))=Cother(s_range_C,IDsum_b_over(2,((length(IDsum_b_over)+1)/2)));
                    [tt1,IDW_o,IDsum_o,tt2,tt3]=RidgeID(Cother,s_range_C(1),s_range_C(end),0,SNR_window,gap_snr,gap,N_adover);
                    for IDW_o_i=0:10
                        if IDsum_o{N_adover+1}(1,IDsum_b_over(2,((length(IDsum_b_over)+1)/2))+IDW_o_i)~=0
                            IDW_o_i=IDsum_b_over(2,((length(IDsum_b_over)+1)/2))+IDW_o_i;
                            break;
                        elseif IDsum_o{N_adover+1}(1,IDsum_b_over(2,((length(IDsum_b_over)+1)/2))-IDW_o_i)~=0
                            IDW_o_i=IDsum_b_over(2,((length(IDsum_b_over)+1)/2))-IDW_o_i;
                            break;
                        end
                    end
                    IDW_ov(s_range_C,IDsum_b_over(2,((length(IDsum_b_over)+1)/2)))=IDW_o{N_adover+1}(s_range_C,IDW_o_i);
                end
            end                
            IDsum_b_over=zeros(2,1); 
        end
    end % End of /for IDsum_b_i=2:length(IDsum_b)/
    [cwt_Area,ratio,ratio_matrix,ParaM_p]=CalArea(C_ov,IDsum_ov,IDW_ov,IDsum_b,D_dist,C50,Area,s_range,s_range_Area,w_range,pmin_number,ratio_matrix_n);    
else
    cwt_Area=0;ratio=0;ratio_matrix=0;ParaM_p=0;
end % End of if overlap==0

function [Lpeak,Rpeak,w_L,w_R,Cother]=isolate(Cinput,id,ID5,ID6,IDW,s_range_C,ratio_matrix_n,s_range_Area,pmin_number)
% /This is new vision split overlap peaks consider min CWT coefficient.

number_p=length(id);
length_c=length(Cinput(1,:));
idd=0;
Lpeak=zeros(5,1);Rpeak=zeros(5,1);
if ratio_matrix_n==1
    eval(['load(''','data result/','Gauss.mat',''');']);
elseif ratio_matrix_n==2
    eval(['load(''','data result/','Cauchy.mat',''');']);
elseif ratio_matrix_n==3
    eval(['load(''','data result/','GauCauchy.mat',''');']);
end

for IDsum_b_oi=1:number_p   
    if sum(IDW(:,id(1,IDsum_b_oi)))>0   %fix 0817
        idd=[idd id(1,IDsum_b_oi)];
    end
end
idd=idd(2:end);number_p=length(idd);
    IDsum_b_0=zeros(length(s_range_C),number_p);
    IDsum_b_L=zeros(length(s_range_C),number_p);IDsum_b_LDepth=zeros(length(s_range_C),number_p);
    IDsum_b_R=zeros(length(s_range_C),number_p);IDsum_b_RDepth=zeros(length(s_range_C),number_p);
    IDsum_b_L0=zeros(length(s_range_C),number_p);
    IDsum_b_R0=zeros(length(s_range_C),number_p);
for s=s_range_C 
    for IDsum_b_oi=1:number_p
        if IDW(s,idd(1,IDsum_b_oi))-ID6(s,idd(1,IDsum_b_oi))==0
            IDsum_b_0(s,IDsum_b_oi)=idd(1,IDsum_b_oi);
        else
            for io=1:100
                if IDW(s,idd(1,IDsum_b_oi))-ID6(s,idd(1,IDsum_b_oi)+io)==0
                    IDsum_b_0(s,IDsum_b_oi)=idd(1,IDsum_b_oi)+io;
                    break;
                elseif IDW(s,idd(1,IDsum_b_oi))-ID6(s,idd(1,IDsum_b_oi)-io)==0
                    IDsum_b_0(s,IDsum_b_oi)=idd(1,IDsum_b_oi)-io;
                    break;
                else
                end                                
            end
        end
        b_R0_check=0;b_L0_check=0;
        for io=1:200
            if IDsum_b_0(s,IDsum_b_oi)+io+1>=length_c
            else
                if Cinput(s,IDsum_b_0(s,IDsum_b_oi)+io)>=0&&Cinput(s,IDsum_b_0(s,IDsum_b_oi)+io+1)<0&&b_R0_check==0
                    IDsum_b_R0(s,IDsum_b_oi)=IDsum_b_0(s,IDsum_b_oi)+io;    % 0817 11:22 fix
                    b_R0_check=1;
                end
                if ID5(s,IDsum_b_0(s,IDsum_b_oi)+io)==-1||io+IDsum_b_0(s,IDsum_b_oi)==length_c-1
                    IDsum_b_R(s,IDsum_b_oi)=IDsum_b_0(s,IDsum_b_oi)+io;
                    IDsum_b_RDepth(s,IDsum_b_oi)=Cinput(s,IDsum_b_0(s,IDsum_b_oi)+io);
                    break;
                end
            end
        end
        for io=1:200
            if IDsum_b_0(s,IDsum_b_oi)-io-1<=0
            else
                if Cinput(s,IDsum_b_0(s,IDsum_b_oi)-io)>=0&&Cinput(s,IDsum_b_0(s,IDsum_b_oi)-io-1)<0&&b_L0_check==0
                    IDsum_b_L0(s,IDsum_b_oi)=IDsum_b_0(s,IDsum_b_oi)-io;  % 0817 11:22 fix
                    b_L0_check=1;
                end
                if ID5(s,IDsum_b_0(s,IDsum_b_oi)-io)==-1||io==IDsum_b_0(s,IDsum_b_oi)-1
                    IDsum_b_L(s,IDsum_b_oi)=IDsum_b_0(s,IDsum_b_oi)-io;
                    IDsum_b_LDepth(s,IDsum_b_oi)=Cinput(s,IDsum_b_0(s,IDsum_b_oi)-io);
                    break;
                end
            end
        end
    end
end
[Lpeak(2,1),Lpeak(3,1),Lpeak(1,1),w_L,h_temp_L,center_temp_L]=funcWH(IDsum_b_0,IDsum_b_L,IDsum_b_L0,IDsum_b_LDepth,1,s_range_Area,D_dist,Zero_dist,C2575,w_range,pmin_number);
[Rpeak(2,1),Rpeak(3,1),Rpeak(1,1),w_R,h_temp_R,center_temp_R]=funcWH(IDsum_b_0,IDsum_b_R,IDsum_b_R0,IDsum_b_RDepth,2,s_range_Area,D_dist,Zero_dist,C2575,w_range,pmin_number);
Lpeak(4,1)=Area(2,w_L);Rpeak(4,1)=Area(2,w_R);
if ratio_matrix_n==1||ratio_matrix_n==2
    Lpeak(5,1)=ratio_matrix_n;Rpeak(5,1)=ratio_matrix_n;
elseif ratio_matrix_n==3
    if w_L<=length(D_dist)/2
        Lpeak(5,1)=1;
    else
        Lpeak(5,1)=2;
    end
    if w_R<=length(D_dist)/2
        Rpeak(5,1)=1;
    else
        Rpeak(5,1)=2;
    end  
end
Cother=Cinput-toC(Gau_C{w_L},Lpeak(1,1),length(Cinput(1,:)),length(Cinput(:,1)))-toC(Gau_C{w_R},Rpeak(1,1),length(Cinput(1,:)),length(Cinput(:,1)));

function [cwt_Area,ratio,ratio_matrix,ParaM_p]=CalArea(C,IDsum,IDW,IDsum_b,D_dist,C50,Area,s_range,s_range_Area,w_range,pmin_number,ratio_matrix_n)

ratio_matrix=zeros(length(s_range),length(IDsum));
ratio=zeros(1,length(IDsum));cwt_Area=zeros(1,length(IDsum));
ParaM_p=zeros(5,length(IDsum_b));
for IDsum_b_i=1:length(IDsum_b)     % IDsum_b_i: Peak i location in IDsum_b list
    s_length=0;s_start=0;
    for s=max(s_range_Area(1),s_range(1)):min(s_range_Area(end),s_range(end))
        if IDW(s,IDsum_b(IDsum_b_i))>0
            s_length=s_length+1;
            if s_start==0
                s_start=s;
            end
        elseif IDW(s,IDsum_b(IDsum_b_i))==0&&s_start~=0
            break;
        end
    end
    if s_start==0
    else
        D_dist_temp=D_dist(s_start:s_start+s_length-1,:);
        [ww,range,output]=pmin(abs(IDW(s_start:s_start+s_length-1,IDsum_b(IDsum_b_i))*ones(1,length(D_dist_temp(1,:)))-D_dist_temp),pmin_number); %fix 0824
        range=range+s_start-1;
        C50_temp=C50(range,ww);
        ratio_matrix(range,IDsum_b(IDsum_b_i))=C(range,IDsum_b(IDsum_b_i))./C50_temp;
        ratio(1,IDsum_b(IDsum_b_i))=mean(ratio_matrix(range,IDsum_b(IDsum_b_i)));
        cwt_Area(1,IDsum_b(IDsum_b_i))=Area(2,ww)*ratio(1,IDsum_b(IDsum_b_i));
        ParaM_p(1:end-1,IDsum_b_i)=[IDsum_b(IDsum_b_i);w_range(ww);ratio(1,IDsum_b(IDsum_b_i));cwt_Area(1,IDsum_b(IDsum_b_i))];
        if ratio_matrix_n==1
             ParaM_p(end,IDsum_b_i)=1;
        elseif ratio_matrix_n==2
             ParaM_p(end,IDsum_b_i)=2;
        else
            if ww>length(D_dist)/2
                ParaM_p(end,IDsum_b_i)=2;
            elseif ww<=length(D_dist)/2
                ParaM_p(end,IDsum_b_i)=1;
            end
        end
    end
end

function [x,ParaM,DrawPeak]=Randmake_peak(Length,peakN,w_bound,h_bound,x,ParaM,DrawPeak,func)
%{ 
ex: [x,ParaM,DrawPeak]=Randmake_peak2(100,10,5,5,0,0,0);
ex: [x,ParaM,DrawPeak]=Randmake_peak2(100,10,5,5,x,ParaM,DrawPeak);
 
Length: data length
PeakN: number of peaks
w_bound: peak gaussian parameter "width" range [-w_bound,w_bound] 
h_bound: peak gaussian parameter "heigth" range [-h_bound,h_bound]
 
x: [3xLength] matrix   [data_no.;data_serious; peak_Area] 
ParaM: [4xpeakN] matrix  [peak_location; peak_width; peak_height; peak_Area]
DrawPeak: [peakNxLength] matrix  [peak1_data; peak2_data; ....;peakN_data]
%}
dd=1;x=0;
while x==0
    while dd==1;
        RandM=rand(3,peakN);
        RandM=sortrows(RandM');RandM=RandM';
        ParaM=[round(RandM(1,:)*Length);RandM(2,:)*w_bound;RandM(3,:)*h_bound;zeros(1,peakN)]; 
        if all(ParaM(1,2:peakN)-ParaM(1,1:peakN-1))~=0
            dd=0;break;
        end
    end
    DrawPeak=zeros(peakN,Length);
    for i=1:peakN
        if exist('x0','var')==0
             x0=[1:Length;zeros(2,Length)];
        end
        [x0,A]=make_peak(x0,ParaM(1,i),ParaM(2,i),ParaM(3,i),0,func);
        if i==1
            ParaM(4,i)=A;
        else
        ParaM(4,i)=A-sum(ParaM(4,1:i-1));
        end
        x0(3,ParaM(1,i))=ParaM(4,i);
        if func==1
            DrawPeak(i,:)=gaussmf(1:Length,[ParaM(2,i) ParaM(1,i)])*ParaM(3,i);
        elseif func==2
            DrawPeak(i,:)=cauchypdf(1:Length,ParaM(1,i),ParaM(2,i))*ParaM(3,i);
        end
    end   
x=x0;
end

function [out,Area]=make_peak(in,center,width,high,pic,func)
% To make a simulated gaussian peak 
% Example: [x0,A]=make_peak(x0,ParaM(1,i),ParaM(2,i),ParaM(3,i),0,func);

%func=1:gaussm function,  func=2: cauchy function

%0 ===== Parameters setup====
if func==1
    out=gaussmf(in(1,:),[width center])*high;
elseif func==2
    out=cauchypdf(in(1,:),center,width)*high;
end
out=out+in(2,:);
out=[in(1,:);out;in(3,:)];
if pic==1
    subplot(2,1,1);plot(in(2,:));
    subplot(2,1,2);plot(out(2,:));
else
end
Area=trapz(out(2,:));

function [snr,noisy_level]=SNR(C,scale,position,window,gap)
%{
ex: [snr,noisy_level]=SNR(C,1:20,0,100);
[ Calculate Signal to Noise level(SNR) ]

snr:  Signal to Noise level(SNR)
C:  Continuous wavelet transform coefficient matrix
scale:  CWT scale
position:  specific m/z position matrix
window:  calculate SNR window  
%}
[a,b]=size(C);
if position==0&scale==0
    scale=1:a;position=1:b;
elseif position==0
    position=1:b;
elseif scale==0
    scale=1:a;
end    
p_length=length(position);
s_length=length(scale);
snr=zeros(s_length,p_length);noisy_level=zeros(s_length,p_length);
for s=scale(1,1):scale(1,s_length);
    for p=position(1,1):position(1,p_length)
        halfwindow=fix(window/2);
        if p<1||p>length(C(s,:))
        elseif p<halfwindow+1
            x0=C(s,1:window);
        elseif p>length(C(s,:))-halfwindow
            x0=C(s,end-window+1:end);
        else
            x0=C(s,p-halfwindow:p+halfwindow-1);
        end
        x1=sort(abs(x0));
        noisy_n=fix(window*gap)+1;
        if noisy_n>=window
            noisy_n=window;
        end
        noisy_level(s,p-position(1,1)+1)=x1(1,noisy_n);
        snr(s,p-position(1,1)+1)=C(s,p)/noisy_level(s,p-position(1,1)+1);
    end
end

function [output1,output2]=Wdist(input,iii,mb,output1)
% [ID{end-1}(s,:),ID{end}(s,i)]=Wdist(ID{end-2}(s,:),i,mb,ID{end-1}(s,:));
   l=iii-1-1;r=iii+1;
   while exist('lend','var')~=1&&l>1
       if input(1,l)==-1;
           output1(1,l)=-1;
           lend=l-1;
       end
       l=l-1;
   end
   while exist('rend','var')~=1&&r<mb
       if input(1,r)==-1;
           output1(1,r)=-1;
           rend=r;
       end
       r=r+1;
   end
   if exist('lend','var')==1&&exist('rend','var')==1
       output2=rend-lend;
   else
       output2=0;
   end
   
function [AB,logic]=match(A,B,adover)
%ex: [AB,logic]=match(fresult,ii,adover)
iitemp=0;
if isempty(A)~=1&&isempty(B)~=1
for a=1:length(A)
    repeat=0;
    for b=1:length(B)
        if abs(A(a)-B(b))<=adover&&repeat==0
            iitemp=[iitemp A(a)];
            repeat=1;            
        end
    end
end
if length(iitemp)==1
    logic=0;
    AB=0;
else
    logic=1;
    AB=iitemp(2:end);
end
else
    logic=0;
    AB=0;
end 

function [w,h,center,ww,h_temp,center_temp]=funcWH(IDsum_b_0,IDsum_b_x,IDsum_b_x0,IDsum_b_xDepth,lr,s_range_Area,D_dist,Zero_dist,C2575,w_range,pmin_number)

[ww,range,output]=pmin(abs(abs(IDsum_b_x(s_range_Area,lr)-IDsum_b_x0(s_range_Area,lr))*ones(1,length(D_dist(1,:)))-(D_dist(s_range_Area,:)-Zero_dist(s_range_Area,:))/2),pmin_number);clear aa;
h_temp=-IDsum_b_xDepth(range,lr)./C2575(range,ww);
w=w_range(1,ww);
center_temp=IDsum_b_0(range,lr);
h=mean(h_temp);center=round(mean(center_temp));

function outC=toC(Gau_C,center,lengC,a)
b=length(Gau_C(1,:));
if center<b/2+1
    outC=[Gau_C(1:a,b/2-center+1:end) zeros(a,lengC-center-b/2)];
elseif center>lengC-b/2-1
    outC=[zeros(a,center-b/2) Gau_C(1:a,1:b/2+lengC-center)];
else
    outC=[zeros(a,center-b/2) Gau_C(1:a,:) zeros(a,lengC-center-b/2)];
end
    
function p= cauchypdf(x, varargin)
a=	0.0;b=	1.0;
if(nargin >= 2)
    a=	varargin{1};
    if(nargin == 3)
        b=			varargin{2};
        b(b <= 0)=	NaN;	% Make NaN of out of range values.
    end
end
if((nargin < 1) || (nargin > 3))
    error('At least one argument, at most three!');
end	
p=	b./(pi*(b.^2 + (x-a).^2));

function [w,range,output]=pmin(input,number)
[sn,n]=size(input);
if sn<number
    number=sn;
end
if number>=sn
    [mina,w]=min(sum(input,1));
    range=1:sn;output=input(range,w);
else
    temp=zeros(2,sn-number+1);
    for s=1:sn-number+1
        [temp(1,s),temp(2,s)]=min(sum(input(s:s+number-1,:),1));
    end 
    [temp_a,s_start]=min(temp(1,:));
    range=s_start:s_start+number-1;
    w=temp(2,s_start);output=input(range,w);
end
    
    