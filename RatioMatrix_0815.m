function RatioMatrix_0815(func)
%ex: RatioMatrix(1);        for inner simulated peaks data (choose1) 
if func==1
    [D_dist,Zero_dist,C50,C2575,Area,Gau_x,Gau_C,s_range,w_range]=RM(func);
    eval(['save(''','data result/','Gauss.mat',''',''D_dist'',''Zero_dist'',''C50'',''C2575'',''Area'',''Gau_x'',''Gau_C'',''s_range'',''w_range'');']);
elseif func==2
    [D_dist,Zero_dist,C50,C2575,Area,Gau_x,Gau_C,s_range,w_range]=RM(func);
    eval(['save(''','data result/','Cauchy.mat',''',''D_dist'',''Zero_dist'',''C50'',''C2575'',''Area'',''Gau_x'',''Gau_C'',''s_range'',''w_range'');']);
elseif func==3
    [D_dist_1,Zero_dist_1,C50_1,C2575_1,Area_1,Gau_x_1,Gau_C_1,s_range,w_range_1]=RM(1);
    [D_dist_2,Zero_dist_2,C50_2,C2575_2,Area_2,Gau_x_2,Gau_C_2,s_range,w_range_2]=RM(2);
    D_dist=[D_dist_1 D_dist_2];Zero_dist=[Zero_dist_1 Zero_dist_2];C50=[C50_1 C50_2];C2575=[C2575_1 C2575_2];Area=[Area_1 Area_2];
    Gau_x=[Gau_x_1;Gau_x_2];Gau_C=[Gau_C_1 Gau_C_2];w_range=[w_range_1 w_range_2];
    eval(['save(''','data result/','GauCauchy.mat',''',''D_dist'',''Zero_dist'',''C50'',''C2575'',''Area'',''Gau_x'',''Gau_C'',''s_range'',''w_range'');']);
end

function [D_dist,Zero_dist,C50,C2575,Area,Gau_x,Gau_C,s_range,w_range]=RM(func)
%% 0 ===== parameter =====
%func=1;                   %func=1 gaussian function  func=2  cachy function  func=3 gaussian+cachy  function
dataLength=1000;
w_range=1:1:50;
s_range=1:1:30;    
mm_window=0;              % mm_window: window size of min/max sliding window
SNR_window=100;           % SNR_window: window size of SNR sliding window
N_adover=2;               % max error adjacent tolerant in over-scales analysis 0,1,2
gap=0;                  % gap: threhold of min/max value of CWT coefficient
gap_snr=0;                % gap_snr: cwt(s,b)/noise_level>gap_snr as meaningful peak

center=round(dataLength/2);
D_dist=zeros(length(s_range),length(w_range));
Zero_dist=zeros(length(s_range),length(w_range));
C50=zeros(length(s_range),length(w_range));
C2575=zeros(length(s_range),length(w_range));
Area=zeros(2,length(w_range));
Gau_x=zeros(length(w_range),dataLength);
Gau_C(1:length(w_range))={zeros(length(s_range),length(w_range))};

%% 11 ===== Main =====
  
for w_no=1:length(w_range)      
    clear ID IDW IDsum snr noisy_level x;
    if dataLength<=SNR_window
        SNR_window=dataLength;
    end
    [x,Area(2,w_no)]=make_peak([1:dataLength;zeros(2,dataLength)],center,w_range(1,w_no),1,0,func);
    Gau_x(w_no,:)=x(2,:);Area(1,w_no)=w_range(1,w_no);
    Gau_C{w_no}=cwt(Gau_x(w_no,:),s_range,'mexh');           
    [ID,IDW,IDsum,snr,noisy_level]=RidgeID3(Gau_C{w_no},s_range(1),s_range(end),mm_window,SNR_window,gap_snr,gap,N_adover);
    c1=0;c2=0;c3=0;c4=0;
    rl0=zeros(length(s_range),2);
    for d=1:round(dataLength/2)-2
        if c1==0
            if sum(IDW{end}(:,center))>0
                D_dist(:,w_no)=IDW{end}(:,center);
                C50(:,w_no)=Gau_C{w_no}(:,center);
                c1=1;
            elseif sum(IDW{end}(:,center-d))>0
                D_dist(:,w_no)=IDW{end}(:,center-d);
                C50(:,w_no)=Gau_C{w_no}(:,center-d);
                c1=1;
            elseif sum(IDW{end}(:,center+d))>0
                D_dist(:,w_no)=IDW{end}(:,center+d);
                C50(:,w_no)=Gau_C{w_no}(:,center+d);
                c1=1;
            end
        end
        if c2<length(s_range)||c3<length(s_range)
            for s=s_range
                if Gau_C{w_no}(s,center-d+1)>=0&&Gau_C{w_no}(s,center-d)<0
                    rl0(s,1)=center-d+1;
                    c2=c2+1;
                end
                if Gau_C{w_no}(s,center+d-1)>=0&&Gau_C{w_no}(s,center+d)<0
                    rl0(s,2)=center+d-1;
                    c3=c3+1;
                end
            end
        end
        if c1>0&&c2>0&&c3>0
            for s=s_range
                if ID{end-1}(s,center-d)==-1
                    C2575(s,w_no)=abs(Gau_C{w_no}(s,center-d));
                    c4=c4+1;
                end
            end            
        end
        if c4==length(s_range)
            break;
        end
    end
    Zero_dist(:,w_no)=rl0(:,2)-rl0(:,1);
end
 
function [out,Area]=make_peak(in,center,width,high,pic,func)
% To make a simulated gaussian peak 
% Example: [out,Area]=make_peak(in,500,20,5);

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

function [ID,IDW,IDsum,snr,noisy_level]=RidgeID3(input,init_scale,end_scale,mm_window,SNR_window,gap_snr,gap,N_adover)
%{
Ex: [ID,IDW,IDsum,snr,noisy_level]=RidgeID3(C,ridgeS_start,ridgeS_end,mm_window,SNR_window,gap_snr,gap,N_adover);
made 20100402 Tzu-Ching Wu

%% < Method parameter> 
mm_window=0;        % mm_window: window size of min/max sliding window
SNR_window=1000;    % SNR_window: window size of SNR sliding window
gap_snr=3;          % gap_snr: cwt(s,b)/noise_level>gap_snr as meaningful peak
gap=0.5;            % threhold of min/max value of CWT coefficient
N_adover=2;         % max error adjacent tolerant in over-scales analysis 0,1,2
%}


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