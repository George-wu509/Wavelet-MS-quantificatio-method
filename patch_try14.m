function patch_try14()

%% 0 ===== Dictionary setup =====
cdn=findstr(cd,'CWT');ccd=cd;
new_fdir0=[ccd(1:cdn+2),'/'];clear cdn ccd;
new_fdir11='data source/';new_fdir12='data result/';new_fdir13='coding sub/';
path(strcat(new_fdir0,new_fdir13),path)
source0=dir([new_fdir0,new_fdir11]);
source_files={source0(3:end,1).name};[a0,a]=size(source_files);clear a0 source0;

%% 1 ===== batch run =====
for d=1:a
   display(['dataset_',int2str(d)]);tic;
   d_be=findstr(source_files{1,d},'_');dd=str2double(source_files{1,d}(d_be+1:end));
   if exist([new_fdir0,new_fdir12,source_files{1,d},'/'],'dir')==0
       mkdir([new_fdir0,new_fdir12,source_files{1,d},'/']);
   end
   source1=dir([new_fdir0,new_fdir11,source_files{1,d},'/RawSpectra/*.txt']);
   source2=dir([new_fdir0,new_fdir11,source_files{1,d},'/truePeaks/*.txt']);
   [a1,a2]=size(source1);
   source_files1={source1(:,1).name};
   source_files2={source2(:,1).name};clear source2 source1 a2      
   for t=1:a1
      t_be=findstr(source_files1{1,t},'y');t_af=findstr(source_files1{1,t},'.');tt=str2double(source_files1{1,t}(t_be+1:t_af-1));
      noisy0=textread([new_fdir0,new_fdir11,source_files{1,d},'/RawSpectra/',source_files1{1,t}],'%s');
      [ar,br]=size(noisy0);clear br;noisy=zeros(2,ar/2-1);
      for i=1:ar/2-1
          noisy(1:2,i)=str2double(noisy0(2*i+1:2*i+2,1));
      end      
      truth0=textread([new_fdir0,new_fdir11,source_files{1,d},'/truePeaks/',source_files2{1,t}],'%s');
      [ar,br]=size(truth0);clear br;truth=zeros(2,ar/2-1);
      for i=1:ar/2-1
          truth(1:2,i)=str2double(truth0(2*i+1:2*i+2,1));
      end 
      [x,RealPeak]=addP(noisy',truth');[xxa,xxb]=size(x);clear xxa
      eval(['save(''',new_fdir0,new_fdir12,source_files{1,d},'/','D',int2str(dd),'x',int2str(tt),'.mat',''',''x'',''RealPeak'');']);
      xx=zeros(1,xxb);
      if log(x(3,:))<=1
          xx(1,:)=0;
      else
          xx(1,:)=log10(x(3,:));
      end
      h1=figure(1);
      subplot(3,1,1);plot(x(2,:));xlabel('m/s');ylabel('Intensity');axis([-inf,inf,-5000,inf]);
      subplot(3,1,2);bar(x(3,:));xlabel('m/s');ylabel('Peak Area');axis([-inf,inf,-inf,inf]);
      h13=subplot(3,1,3);bar(xx(1,:));xlabel('m/s');ylabel('Log10 Area');axis([-inf,inf,-inf,inf]);%set(h13,'YScale','log');
      eval(['saveas(h1,''',new_fdir0,new_fdir12,source_files{1,d},'/','D',int2str(dd),'x',int2str(tt),'.fig',''');']);
   end
   toc;
end
