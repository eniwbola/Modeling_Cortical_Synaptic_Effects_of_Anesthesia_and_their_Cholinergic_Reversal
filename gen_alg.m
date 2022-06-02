clear
format long

jit_bool=1;
jitter_num=1;
  
try
   exc_ccg_cell_jitter ;
catch
exc_ccg_cell_jitter =cell(5,6,5);
end
    
fig_iter=0;
height_diff=[];

dir_locs_pre=pwd;
dir_locs_pre=strcat(dir_locs_pre,'\')
%dir_locs_pre=strcat(dir_locs_pre,'/')
%%
same_list_rem=0;

%%
dir_load_order=cell(1,20);%cell for each of parameter
dir_locs_pre_gks=strcat(dir_locs_pre,'gks*'); %make string with gks extenstion
b=cell(length(dir_locs_pre_gks),20);

dir_locs=dir(dir_locs_pre_gks); % grab all driectors with gks prefix        
              

 for ggg=1:length(dir_locs) % delete the directories with less than 1000 spike time
               pre_pre_dir_string=dir_locs_pre;
          pre_dir_string=strcat(dir_locs_pre,dir_locs(ggg).name);
          try
          rast_exc_pre=load(strcat(pre_dir_string ,'/raster_datexc.txt')) ;
          rast_inh_pre=load(strcat(pre_dir_string ,'/raster_datinh.txt')) ;
        
         catch
            rmdir(pre_dir_string,'s')
             continue
          end
         
       
           
        if (sum(size(rast_inh_pre))<1) || (sum(size(rast_exc_pre))<1)
           rmdir(pre_dir_string,'s')
                continue
        end

        if  (max(rast_inh_pre(:,2))<1000) || (max(rast_exc_pre(:,2))<1000)
            rmdir(pre_dir_string,'s')
           continue
       end
           
           
           p
           
          
 end

dir_locs=dir(dir_locs_pre_gks);

 for ggg=1:length(dir_locs)
               pre_pre_dir_string=dir_locs_pre;
          pre_dir_string=strcat(dir_locs_pre,dir_locs(ggg).name);
         
          parameter_cur = regexp(dir_locs(ggg).name,'\d*\.\d*|(?<=run_)\d*|(?<=mult_)\d*|(?<=_)\d*','Match');
          % grab parameter numbers
          b(ggg)={dir_locs(ggg).name};
          
          
 end
 dir_locs=dir(dir_locs_pre_gks); %be edit 3/7/2022
 parameter_cur_2 = regexp(b,'\d*\.\d*|(?<=run_)\d*|(?<=wi_mult)\d*|(?<=_)\d*','Match');
              
              parameter_cur_3 = {cat(1,parameter_cur_2{:})}  ; 
              hj=parameter_cur_3{:};
              
              unique(parameter_cur_3{:});
                
              for mm=1:size(hj,2)
                 un_param{mm}= unique(hj(:,mm));
                  
              end
              [row col]=find(strcmp(parameter_cur_3{:},'0.0'));
              
             % return
            close all
               find(strcmp(un_param{1},'0.0'));
             %length(un_param{1})
           raster_sub_plot=1;
           ccg_sub_plot=1;
           stim_sub_plot=1;
           zoom_raster_save=0;
           trace_sub_plot_save=0;
           zoom_raster_save=0;
           raster_save=0;
           
           ccg_exc_mag=zeros(1,length(dir_locs));
           ccg_inh_mag=zeros(1,length(dir_locs));
        
           exc_neuron_count=400;
             tic;
             ccg_com_prob=1.0;
             
            t_stim=6000;
            fig_final=1;
            count_bin_size=1;
            mpc_vec=[];
            stim_complexity=0;
            no_stim_complexity=0;
            compute_mpc=1;
            compute_int_comp=1;
            int_com_calc_tony=1;
            fig_ind_order=[];
            time_bin_pat=[0:1:6000];%max(times_full)]% 1 s

            time_bin_pat=[0:2:6000];%max(times_full)] %2 s time bin
      time_bin_pat=[0:1:6000];%max(times_full)] %2 s time bin
      
      
      ideal_ccg_mat=[];
      ideal_inh_ccg_mat=[];
      ideal_raster_cell=cell(1,2);
      
       save_ccg_fig=1;
      save_sing_ccg=1;
      
      
      
     
      auto_initialize=0;
      
     jitter_half_width=20; 
     ref_exc_neuron_total=40;
     ref_inh_neuron_total=40;

      fig_rows_vec=[];
      fig_cols_vec=[];
      anes_vec=[];
      gks_vec=[];
      isi_vec=zeros(length(dir_locs),1)';
      med_freq_vec=zeros(length(dir_locs),1)';
      CI_vec=zeros(length(dir_locs),1)';
      EI_vec=zeros(length(dir_locs),1)';
      mpc_vec=zeros(length(dir_locs),1)';
      mpc_std_vec=zeros(length(dir_locs),1)';
      
      param_mat=[];
        inh_thresh=1.12;%2.00
        
        dir_tot=length(dir_locs);
        total_param_nums=30;
        base_param_nums=20;
      dir_locs_pre=strcat(dir_locs_pre,'/')
      % mpc % int % complexity
% % % %       
% % % %       for qqq=1:length(dir_locs)
% % % %          pre_dir_str=strcat(dir_locs_pre,dir_locs(qqq).name) ;
% % % %          try
% % % %           rast_exc_pre=load(strcat(pre_dir_string ,'/raster_datexc.txt')) ;
% % % %         rast_inh_pre=load(strcat(pre_dir_string ,'/raster_datinh.txt')) ;
% % % %         
% % % %          catch
% % % %             
% % % %              continue
% % % %          end
% % % %         if  max(rast_inh_pre(:,2))<1000
% % % %             
% % % %         end
% % % %             
% % % %       end
      
      
for qqq=1  : length(dir_locs) %=["-.15"]
                                 
    rast_exc=[];
    rast_inh=[];
    stim_dat=[];
    filt_num=10;
                                        
    fig_iter=fig_iter+1;

    pre_pre_dir_string=dir_locs_pre;
    pre_dir_string=strcat(dir_locs_pre,dir_locs(qqq).name);%grab the current director
    parameter_cur = regexp(dir_locs(qqq).name,'\d*\.\d*|(?<=run_)\d*|(?<=wi_mult)\d*|(?<=_)\d*','Match');
    item_str=dir_locs(qqq).name;% grab the current directory name
    
    cur_run=regexp(pre_dir_string,'(?<=run_)\d*','Match');
    if iscell(cur_run)      
        cur_run='1';% cur_run{2};
    end
                  
    if iscell(cur_run)
       cur_run='1';%cur_run{2};
    end
    
    if( cur_run=='1')
                 
        run_num=0;
        
        for nnn=1%:3 %1:1
    
            pa = regexp(pre_dir_string,'(?<=run_)\d*','Match');
            if iscell(pa)
                pa='1';%pa{2}
            end
    
        pre_dir_string= strrep(pre_dir_string,strcat('run_',pa),strcat('run_',num2str(nnn)));
        
        output_dir_string=strcat(pre_dir_string,'/fig'); %strcat('C:\Users\eniwbola\Desktop\Research\Hudetz\for_git_new\spatially_connected\Michal_Edit\no_stim_scan\fig3');%stim_',stim_range(ooo),'\AE_',AE_range(iii),'\AI_',AI_range(jjj),'\TauS_I_',TauS_I_range(mmm),'\WI_',WI_range(lll),'\igks_',gks_range(hhh),'_RUN_',Run_range(nnn));
  
        output_dir_string_2=strcat(pre_pre_dir_string,'/fig');
        
%------------------------------------------------------------------------------------------------------

        label_index_0=11;%14 %pee
        label_index_1=10;%IE   
        label_index_2=1;%10%22%1;%10; %16;%gks %noise_freq
        label_index_3=14;%8%14%8%14%23%14%6%14%8%14%6%14%6%14%6 %8 % 14;%8;%6;%6;%14; %6;%6%1; %8; %14 %stim_strength %wiE
        output_dir_org=strcat(pre_pre_dir_string,'/wei_',parameter_cur{label_index_0});

        parameter_name='|W:';
        title_add=strcat('gks:',parameter_cur(label_index_2),parameter_name,  un_param{label_index_3}{find(strcmp(un_param{label_index_3}, parameter_cur{label_index_3}))} ); %,'|IE:',parameter_cur(label_index_1),'|pee:',parameter_cur(label_index_0) );

        bd=strcat('we',parameter_cur{label_index_0});
        if(iscell(bd))
            bd=bd{1};
        end

save_add_ccg=strcat('we',parameter_cur{label_index_0});
save_add_raster=strcat('we',parameter_cur{label_index_0});



%end of paramter choice
%------------------------------------------------------------------------%
if(iscell(output_dir_org))
output_dir_org=output_dir_org{1};
end

fig_row_tot=length(un_param{label_index_2});
fig_col_tot=length(un_param{label_index_3});

fig_rc_num= find(strcmp(un_param{label_index_0},parameter_cur{label_index_0}))*10 +find(strcmp(un_param{label_index_1},parameter_cur{label_index_1})) ;
fig_rc_tot= length(un_param{label_index_0})*length(un_param{label_index_1}) ;

fig_row=find(strcmp(un_param{label_index_2},parameter_cur{label_index_2}));
fig_col=find(strcmp(un_param{label_index_3},parameter_cur{label_index_3}));
anes_vec=[anes_vec, str2num(parameter_cur{label_index_3}) ];
gks_vec=[gks_vec , str2num(parameter_cur{label_index_2})];
fig_cols_vec=[fig_cols_vec, fig_col];
fig_rows_vec=[fig_rows_vec, fig_row];
fig_ind_num=(fig_row-1)*fig_col_tot + fig_col;
fig_ind_order=[fig_ind_order; fig_ind_num];

fig_cont_pre_0=find(strcmp(un_param{label_index_0},parameter_cur{label_index_0}));
fig_cont_pre_1=find(strcmp(un_param{label_index_1},parameter_cur{label_index_1}));

fig_cont=(fig_cont_pre_0-1)*length(un_param{label_index_1})+fig_cont_pre_1;





if((qqq==1) )

    
    isi_mat=zeros(fig_row_tot,fig_col_tot,fig_rc_tot);

ccg_isi=zeros(fig_row_tot,fig_col_tot,fig_rc_tot);


mpc_mat=zeros(fig_row_tot,fig_col_tot,fig_rc_tot); %[ mpc_vec mean_mpc];
mpc_std_mat=zeros(fig_row_tot,fig_col_tot,fig_rc_tot);
integration_mat=zeros(fig_row_tot,fig_col_tot,fig_rc_tot);
complexity_mat=zeros(fig_row_tot,fig_col_tot,fig_rc_tot);
CI_mat=zeros(fig_row_tot,fig_col_tot,fig_rc_tot);
EI_mat=zeros(fig_row_tot,fig_col_tot,fig_rc_tot);
integration_mat_no_zero=zeros(fig_row_tot,fig_col_tot,fig_rc_tot);
complexity_mat_no_zero=zeros(fig_row_tot,fig_col_tot,fig_rc_tot);

integration_mat_pre=zeros(fig_row_tot,fig_col_tot,fig_rc_tot);
complexity_mat_pre=zeros(fig_row_tot,fig_col_tot,fig_rc_tot);

integration_mat_post=zeros(fig_row_tot,fig_col_tot,fig_rc_tot);
complexity_mat_post=zeros(fig_row_tot,fig_col_tot,fig_rc_tot);

raster_neuron_cell=cell(fig_row_tot,fig_col_tot);
raster_time_cell=cell(fig_row_tot,fig_col_tot);
raster_inh_neuron_cell=cell(fig_row_tot,fig_col_tot);
raster_inh_time_cell=cell(fig_row_tot,fig_col_tot);




%-----------------------------------------

end

g=[];
for i =1:23
    g=[g str2num(parameter_cur{i})];
end

param_mat=[param_mat;g];




title_string=title_add; %strcat('gks-',parameter_cur{1},'|net-',parameter_cur{2},'|WE-',parameter_cur{6},'|WI',parameter_cur{4},'|conp-',parameter_cur{5},'-N_str-',parameter_cur{8},'-N_freq-',parameter_cur{8});                             

suf_string=strcat('gks_',parameter_cur{1},'_net_',parameter_cur{2},'_WE_',parameter_cur{6},'_WI_',parameter_cur{4},'_conp_',parameter_cur{5},'_N_str_',parameter_cur{8},'_N_freq_',parameter_cur{8},'.png');                             

dir_suf=suf_string;
title_suf=title_string;
mkdir(output_dir_string);
mkdir(output_dir_string_2);
mkdir(output_dir_org);


movie_dir_string=output_dir_string;
fig_dir_string=output_dir_string;


                                       end
                                       
                                   end     
                                       



dir_load_order{qqq}=pre_dir_string;


%-------------------------------begin_int_comp_mpc----------------------%


%%---------- beginning pattern-----------------------------
try
rast_exc_pre=load(strcat(pre_dir_string ,'/raster_datexc.txt')) ;
rast_inh_pre=load(strcat(pre_dir_string ,'/raster_datinh.txt')) ;
catch
 
   CI_vec(qqq)=1000;
   EI_vec(qqq)=1000;
   mpc_vec(qqq)=1000;
   freq_vec(qqq)=1000;
continue

end
  exc_neurons_rast=rast_exc_pre(:,1);
   exc_times_rast=rast_exc_pre(:,2);
   inh_neurons_rast=rast_inh_pre(:,1);
   inh_times_rast=rast_inh_pre(:,2);

neurons_full=exc_neurons_rast;
times_full=exc_times_rast;


isi_vec(qqq)=max(exc_times_rast)/(length(exc_times_rast)/max(exc_neurons_rast)) ;
new_freq=[];
for lk=1:max(exc_neurons_rast)
    sing_times=exc_times_rast(find(exc_neurons_rast==lk));
    isi=max(sing_times)/length(sing_times);
    
new_freq=[ new_freq (1./(isi./1000))];
end
med_freq_vec(qqq)=median(new_freq);
%% setting up patterns

  %% so this the else if we dont have stim complexity

  
m_t=median(times_full);
m_t_p=m_t+2000;
m_t_m=m_t-2000;
  
 m_t_m=600;
  neurons_part=neurons_full.* (  (times_full>(m_t_m) ) ) .*( (times_full<(m_t_p) )     )    ;
times_part=times_full.*(  (times_full>(m_t_m) ) ) .*( (times_full<(m_t_p) )     )  ;  


  neurons_part=neurons_part(neurons_part~=0);
 times_part=times_part(times_part~=0);




neuron_types=unique(neurons_full);
patterns_all=[];
patterns_all_int=[];% this is saving the patterns as a binary matrix instread of as strings so I can do more with it


pat_neurons=randi([1 exc_neuron_count],1,60);

%ismember(neurons_full,pat_neurons)
new_neurons=neurons_part(ismember(neurons_part,pat_neurons));
new_times=times_part(ismember(neurons_part,pat_neurons));
bin_size=1;
cell_it=0;
for hh=1:(round(length(time_bin_pat)/bin_size)-1)
    cell_it=cell_it+1;
    neurons_spiked=new_neurons(find(  (new_times>time_bin_pat(hh)).*(new_times<time_bin_pat(hh+bin_size))));
    new_neurons_un=unique(new_neurons);
    neuron_types=ismember(new_neurons_un,neurons_spiked); % vector with 1's in the bplace of the spiked neurons
    patterns_all_int=[patterns_all_int;neuron_types'];
    pattern_vect=num2str(neuron_types);
    pattern_pre=pattern_vect(find(~isspace(pattern_vect)));
    pattern=convertCharsToStrings(pattern_pre);
    patterns_all=[patterns_all;pattern];
    patterns_2{cell_it}=neurons_spiked;
  %  if .>(hh/length(time_bin_pat)) %(i>length(time_bin_pat)/10) 
    
   % end
   
    hh/length(time_bin_pat);
end



%[pat_un, BB, frequen]=unique(patterns_all)

pat_un=unique(patterns_all);
pat_freq=[];
pat_spikes=[];
for v=1:length(pat_un)
pat_freq=[pat_freq sum(patterns_all==pat_un(v))];
bv=num2str(pat_un(v));
bv=bv((bv)~='0');
pat_spikes=[pat_spikes length(bv)];
end

%unique(patterns_all)
                       

%%------------end patterns-----------------------------%%
num_neurons_mpc=60;
%num_neurons_mpc=400;
rast_exc=rast_exc_pre;
rast_inh=rast_inh_pre;

rast_neu=rast_exc(:,1);
rast_ti=rast_exc(:,2);
times_exc_10=rast_ti(rast_neu<(num_neurons_mpc+1));
neurons_exc_10=rast_neu(rast_neu<(num_neurons_mpc+1));

neurons_exc_10=neurons_exc_10(times_exc_10>1000);
times_exc_10=times_exc_10(times_exc_10>1000);

rast_10=[times_exc_10(:),neurons_exc_10(:)];


if compute_int_comp


if(int_com_calc_tony)
com_string="complexity"
qqq
 try
     [TC,EI,CI, ALLEn, psp]=unitentroJAV3nozero(patterns_all_int);
 catch
     TC=0;
     EI=1000;
     CI=1000;
     ALLEn=0; 
     psp=0;
     %continue
 end
     %[complexity_val]=complexity_fun(patterns_all_int);
   % EI_mat(fig_row,fig_col,fig_cont)=EI;
  EI_vec(qqq)=EI;
   % CI_mat(fig_row,fig_col,fig_cont)= CI;
  CI_vec(qqq)=CI;
   
    
end

integration_mat_no_zero(fig_row,fig_col,fig_cont)=EI;
complexity_mat_no_zero(fig_row,fig_col,fig_cont)=CI;


                             end



com_string="mpc"
qqq
if compute_mpc
   try    
   [mean_mpc,mpc_cellpairs]=mpc_network_class(num_neurons_mpc,rast_10);
        
    mpc_vec(qqq)=mean_mpc;
    mpc_std_vec(qqq)=std(mpc_cellpairs);
   catch
       
        
    mpc_vec(qqq)=1000;
    mpc_std_vec(qqq)=2;
   end
    
end


%--------------------------------- end mpc_int_cont----------------------%

                             end
                                 runtime=   toc
             
%%
file_name='iter_text.txt' ;     
if isfile(file_name)
    fid=fopen(file_name,'r');
    cur_iter=fscanf(fid,'%d')  ;
    cur_iter=cur_iter+1;
    fclose(fid);
    fid=fopen(file_name,'w+');
    fprintf(fid,'%d',cur_iter);
    fclose(fid);
else
   fid=fopen(file_name,'w+');
    cur_iter=1;
    fprintf(fid,'%d',cur_iter);
    fclose(fid);
end
% be edit so here i need to say only if not eaquato to 30
freq_vec=1./(isi_vec./1000);
id_vec=[1:1:length(CI_vec)];
it_vec=cur_iter*ones(1,length(CI_vec));%[1:1:length(CI_vec)];
metric_vec=[it_vec',id_vec',CI_vec',EI_vec',mpc_vec',freq_vec'];
A=metric_vec;

 %dir_tot=length(dir_locs);
 %       total_param_nums=30;
 %       base_param_nums=20;
 if dir_tot>20%total_param_nums
    fid = fopen('Mymatrix.txt','wt');
    for ii = 1:size(A,1)
        fprintf(fid,'%g\t',A(ii,:));
        fprintf(fid,'\n');
    end
    fclose(fid);
          
end


A=metric_vec;
fid = fopen('metrics_current.txt','wt');
for ii = 1:size(A,1)
    fprintf(fid,'%g\t',A(ii,:));
    fprintf(fid,'\n');
end
fclose(fid);


 
  if dir_tot>20 %total_param_nums
A=metric_vec;
fid = fopen('metrics_total.txt','a');
for ii = 1:size(A,1)
    fprintf(fid,'%g\t',A(ii,:));
    fprintf(fid,'\n');
end
fclose(fid);
  end


  
 
  

%% value

%standard_vec=[4.5 .6 .13 15];
standard_vec=[3.5 .62 .05 6.5];

%ratio_vec= metric_vec(:,2:end)./standard_vec

ratio_vec=((metric_vec(:,3:end)-standard_vec) ./standard_vec).^2;
cost_vec=sum(ratio_vec,2);%sqrt(sum(ratio_vec,2))
it_vec=cur_iter*ones(1,length(CI_vec));%[1:1:length(CI_vec)];
id_vec=[1:1:length(CI_vec)];
cost_mat=[it_vec',id_vec',ratio_vec,cost_vec]

A=cost_mat;

 
fid = fopen('cost_current.txt','wt');
for ii = 1:size(A,1)
    fprintf(fid,'%g\t',A(ii,:));
    fprintf(fid,'\n');
end
fclose(fid);



  if dir_tot>20%total_param_nums
fid = fopen('cost_total.txt','a');
for ii = 1:size(A,1)
    fprintf(fid,'%g\t',A(ii,:));
    fprintf(fid,'\n');
end
fclose(fid)
  end
%%
id_vec=[1:1:length(CI_vec)];
param_mat_2=[it_vec',id_vec',param_mat]


A=param_mat_2;

fid = fopen('param_mat_current.txt','wt');
for ii = 1:size(A,1)
    fprintf(fid,'%g\t',A(ii,:));
    fprintf(fid,'\n');
end
fclose(fid);


  if dir_tot>20%=total_param_nums
A=param_mat_2;
fid = fopen('param_mat_total.txt','a');
for ii = 1:size(A,1)
    fprintf(fid,'%g\t',A(ii,:));
    fprintf(fid,'\n');
end
fclose(fid);
  end
 
  
% param_str=["gks","pee","pii","pie","pei","wee","wii","wie","wei","CIe","CIi","stim","w_mult","nmda_wee","run","nmda_wei","log_mu","log_std","n_freq","n_sten","n_dur","net","run_tim"]
% 'gks_1.2_pee_0.1_pii_0.1_pie_0.1_pei_0.1_wee_1.2_wii_1.4_wie_1.4_wei_1.2_CIe_0.0_CIi_0.0_stim_0.0_w_mult_1_nmda_strength_wee_0.37_run_1_nmda_strength_wei_0.37_log_mu_-20.0_log_std_9.4_noise_freq_0.1_noise_strength_4.0_noise_dur_2.0_net_10_ruin_time_20000'

%%
param_str=["gks","pee","pii","pie","pei","wee","wii","wie","wei","CIe","CIi","stim","w_mult","nmda_wee","run","nmda_wei","log_mu","log_std","n_freq","n_sten","n_dur","net","run_tim"]
% 'gks_1.2_pee_0.1_pii_0.1_pie_0.1_pei_0.1_wee_1.2_wii_1.4_wie_1.4_wei_1.2_CIe_0.0_CIi_0.0_stim_0.0_w_mult_1_nmda_strength_wee_0.37_run_1_nmda_strength_wei_0.37_log_mu_-20.0_log_std_9.4_noise_freq_0.1_noise_strength_4.0_noise_dur_2.0_net_10_ruin_time_20000'
%param_mat(22,[2 3 4 5 1 6 7 8 9 10 11 15 12 13 14 16 17 18 19 20 21 22 23])
%loc_vec=[2 3 4 5  6 7 8 9  10 11 1  12 15  13 14  16 17 18 19 20 21 22 23]
loc_vec=[2 3 4 5  6 7 8 9  10 11 1  12 13 15  14  16 17 18 19 20 21 22 23]

% for ll=1:length(loc_vec)
% strcat(param_str(loc_vec(ll)),num2str(param_mat(22,loc_vec(ll)) ))
% 
% end

%% parmater choice

%have twenty sim bottom ten with changes and then cmpare to rest and
%designate new twenty

% so have thirty, remove lowest ten and resimulate mid ten and then get
% thirrty again 

% so in practice first sim thirty, 
% decide lowest ten and save params from mid ten
% put lowest ten in box and and randomly purmut paramters betweeen mid and
% highest ten with both an random averaging and  and a random variable
% that is is between is the distance

if dir_tot>=30
      dir_diff=0;
else
    
  dir_diff=30-dir_tot;
end

[val ind]=sort(cost_mat(:,3));

if dir_tot>=30
    
ind_low=ind(1:(10-dir_diff));ind_mid=ind((11-dir_diff):(20-dir_diff));
ind_high=ind((21-dir_diff):(30-dir_diff)); ind_last=ind(30-dir_diff+1:end);
elseif dir_tot>20 && (dir_tot<=30)
ind_low=ind(1:(10-dir_diff));ind_mid=ind((11-dir_diff):(20-dir_diff));
ind_high=ind((21-dir_diff):(30-dir_diff)); ind_last=ind(30-dir_diff+1:end);   
elseif dir_tot>=10  && (dir_tot<=20)
ind_low=ind(1:10);ind_mid=ind(11:end);
ind_high=[]; ind_last=[];    
elseif (dir_tot>=0) && (dir_tot<10)
ind_low=ind(1:end);ind_mid=[];
ind_high=[]; ind_last=[];

end


param_num=round(size(param_mat,1));    

half_param=floor(param_num/2);
ind_inter_low=ind(1:(half_param)); %ind(1:(10-dir_diff))
ind_inter_mid=ind((half_param)+1:2*half_param);  %ind(1:(10-dir_diff))



dir_locs_high={};
dir_locs_mid={};
dir_locs_low={};
for m =1:length(ind_high)
dir_locs_high{m}=dir_locs(ind_high(m)).name;
end

for m =1:length(ind_mid)
dir_locs_mid{m}=dir_locs(ind_mid(m)).name;
end

for m =1:length(ind_low)
dir_locs_low{m}=dir_locs(ind_low(m)).name;
end

% ind_last=ind(30:end);

for m =1:length(ind_last)
dir_locs_last{m}=dir_locs(ind_last(m)).name;
end

% the

% mkdir('rejects')
% for m =1:length(ind_last)
% movefile(dir_locs_last{m}, 'rejects')
% end


%if dir_tot>=total_param_nums
%if dir_tot>20%total_param_nums% put 20 just to trobul shouot 
if dir_tot>=total_param_nums
mkdir('rejects');
try
for m =1:length(ind_last)
try
movefile(dir_locs_last{m}, 'rejects');
catch
%delete(dir_locs_high{m}) 
continue
end

end
catch
   h="hey" ;
end

for m =1:length(ind_high)
try
movefile(dir_locs_high{m}, 'rejects');
catch
%delete(dir_locs_high{m}) 
continue
end

end

end



%% choose new parameters

% % % mid_params=param_mat(ind_mid,2:end);
% % % low_params=param_mat(ind_low,2:end);
% % % 
% % % low_mid_dif=(low_params-mid_params);
% % % 
% % % new_params=mid_params+low_mid_dif/2;
% % %   
% % % new_params=new_params+ 1.5*(rand-.5)*low_mid_dif;
% % 
% % %if dir_tot>=total_param_nums
% % if dir_tot>20
% % mid_params=param_mat(ind_mid,loc_vec);
% % low_params=param_mat(ind_low,loc_vec);
% % 
% % low_mid_dif=(low_params-mid_params);
% % else
% %     
% % %     ind_inter_mid
% % %     ind_inter_low
% %    mid_params=param_mat(ind_inter_mid,loc_vec);
% % low_params=param_mat(ind_inter_low,loc_vec);
% % 
% % low_mid_dif=(low_params-mid_params); 
% %     
% % % param_num=round(size(param_mat,1));    
% % % p = randperm(param_num);
% % % half_param=floor(param_num/2);
% % % mid_params=param_mat((1:half_param),loc_vec);
% % % low_params=param_mat((half_param+1):(half_param*2),loc_vec);
% % % 
% % % low_mid_dif=(low_params-mid_params);   
% %     
% %     
% % end
% % 
% % 
% % %new_params=mid_params+low_mid_dif/2;
% %   
% % mod_location=[4,5,6,7,15]; %which parameters I want to modify
% % 
% % new_params(:,mod_location)=new_params(:,mod_location)+ 1.5*(rand-.5)*low_mid_dif(:,mod_location);
%% So i need to do the above but ind random 

% mid_params=param_mat(ind_mid,2:end);
% low_params=param_mat(ind_low,2:end);
% 
% low_mid_dif=(low_params-mid_params);
% 
% new_params=mid_params+low_mid_dif/2;
%   
% new_params=new_params+ 1.5*(rand-.5)*low_mid_dif;

%if dir_tot>=total_param_nums
perm_5_1=nchoosek(ind_mid,min(5,length(ind_mid)));
perm_5_2=nchoosek(ind_low,min(5,length(ind_low)));
my_rand_perm_5_1=randi(size(perm_5_1,1));
my_rand_perm_5_2=randi(size(perm_5_2,1));
%low_inds=[ind_mid(:);ind_low(:)];
low_inds=[perm_5_1(my_rand_perm_5_1,:)';perm_5_2(my_rand_perm_5_2,:)'];
% % low_inds_inter=[ind_inter_mid(:);ind_inter_low(:)];
% % new_param_temp_1=param_mat(low_inds,loc_vec);
% % new_param_temp_2=param_mat(low_inds_inter,loc_vec);
%if dir_tot>20
if dir_tot>=total_param_nums
mid_params=param_mat(ind_mid,loc_vec);
low_params=param_mat(ind_low,loc_vec);

low_mid_dif=(low_params-mid_params);
else
    
%     ind_inter_mid
%     ind_inter_low
mid_params=param_mat(ind_inter_mid,loc_vec);
low_params=param_mat(ind_inter_low,loc_vec);

low_mid_dif=(low_params-mid_params); 
    
% param_num=round(size(param_mat,1));    
% p = randperm(param_num);
% half_param=floor(param_num/2);
% mid_params=param_mat((1:half_param),loc_vec);
% low_params=param_mat((half_param+1):(half_param*2),loc_vec);
% 
% low_mid_dif=(low_params-mid_params);   
    
    
end


%new_params=mid_params+low_mid_dif/2;
  

%mod_location=[4,5,6,7,15]; %which parameters I want to modify
mod_location=[5,6,15,11]; %which parameters I want to modify

%new_params(:,mod_location)=new_params(:,mod_location)+ 1.5*(rand-.5)*low_mid_dif(:,mod_location);


mut_rate=.8;

%new_params=param_mat(sing_param(1:10))
loc_vec_sim=[2 3 4 5  6 7 8 9  10 11 1  12 13 15  14  16 17 18 19 20 21 22 23]

new_params=[];
old_param_mat_pre=[];
old_param_mat_post=[];
for i=1:9
    F=2*rand;
    % pick x and a,b,c
    perm_inds=perms(low_inds);
    rand_loc=randi(length(perm_inds));  
    sing_perm=perm_inds(rand_loc,:);
    %my_rand_i=randi(len(low_inds));
    
    %x=param_mat(sing_perm(1),:);%a=param_mat(sing_perm(2),:);%b=param_mat(sing_perm(3),:);%c=param_mat(sing_perm(4),:);
    x=param_mat(sing_perm(1),loc_vec);
    
    a=param_mat(sing_perm(2),loc_vec);
    b=param_mat(sing_perm(3),loc_vec);
    c=param_mat(sing_perm(4),loc_vec);
    
    
    for i=mod_location
        my_rand=rand;
        if my_rand<mut_rate
            x(i)=a(i)+F*(b(i)-c(i));
        else
            hey="hey";
        end
        if (i==6) || (i==15)
           x(i+1)=x(i) ;
        end
    end
    
    old_param_mat_pre=[old_param_mat_pre;param_mat(sing_perm(1),:)];
    old_param_mat_post=[old_param_mat_post;param_mat(sing_perm(1),loc_vec)];
    new_params=[new_params;x];
end

old_param_mat_pre
old_param_mat_post
new_params
%%



A=new_params;
fid = fopen('new_param.txt','wt');
for ii = 1:size(A,1)    
    fprintf(fid,'%g\t',A(ii,:));
    fprintf(fid,'\n');
end
fclose(fid);






 if dir_diff>0
 try
 A=new_params(1:dir_diff,:)
 catch
      A=new_params;
 end
%fid = fopen('remaining_param.txt','wt');
fid = fopen('new_param.txt','wt');
for ii = 1:size(A,1)    
    fprintf(fid,'%g\t',A(ii,:));
    fprintf(fid,'\n');
end
fclose(fid);
 
 end
%% param selection

%newdir=${cwd}"/gks_"${gks}"_pee_"${pee}"_pii_"${pii}"_pie_"${pie}"_pei_"${pei}"_wee_"${wee_mult}"_wii_"${wii_mult}"_wie_"${wie_mult}"_wei_"${wei_mult}"_CIe_"${CIe}"_CIi_"${CIi}"_stim_"${stim}"_w_mult_"${w_mult}"_nmda_strength_wee_"${wee_nmda_mult}"_run_"${run}"_nmda_strength_wei_"${wei_nmda_mult}"_log_mu_"${log_norm_mu}"_log_std_"${log_norm_std}"_noise_freq_"${noise_freq}"_noise_strength_"${noise_strength}"_noise_dur_"${noise_dur}"_net_"${net}"_ruin_time_"${run_time}"_frac_"${stim_frac}

%% mpc classification

              
   for pp=1:1:length(fig_ind_order) 
    % isi_mat_t=isi_mat';
 %ccg_exc_mat(fig_row,fig_col,fig_cont)=mean(exc_strength_no_nan);%sum(exc_dif_vec);
mpc_mat(fig_rows_vec(pp),fig_cols_vec(pp))=mpc_vec(pp);%sum(exc_dif_vec)
CI_mat(fig_rows_vec(pp),fig_cols_vec(pp)) = CI_vec(pp);%sum(exc_dif_vec)
EI_mat(fig_rows_vec(pp),fig_cols_vec(pp)) =EI_vec(pp);
freq_map(fig_rows_vec(pp),fig_cols_vec(pp))=1./(isi_vec(pp)./1000);
   end
  
  
  
  r=1;
  
  
  return  
  
  
%% start of new connectivity
                  
%-------------------------------------------------- classification  end --------------------------%                              
    %100 jitters 100 connections =1 hour
    

%% unitentroJAV
function [TC, EI, CI, ALLEn, psp]=unitentroJAV3nozero(qf)
%----- unitentro --------
%unitary analysis and entropies
%called from unitary_main

%data should be in matrix qf
%data are 0 (no spike) and 1 (spike)
%columns: channels, rows: time points (e.g. 1 ms bins)

%modified for EI to omit all-zero spike patterns, AGH 4-9-14
%note that this conflicts with CI2 calculation
%added TC2 and CI2 from randomized or shuffled data, AGH 12-10-14 

%---------unique spike patterns--------------------------------------------   
    [R,C]=size(qf);
    qfs=num2str(qf,2); M=cellstr(qfs); %cell array of spike patterns as strings
    [~, ind1]=unique(sort(M));
    count1=diff([0;ind1]); %counts of unique patterns
    count1(1)=[]; %remove no spike pattern ********************************
    psp=count1/sum(count1); %probability of unique spike patterns
    npat=length(count1)-1; %number of unique patterns
   
%---------Unitary events---------------------------------------------------
%     p=sum(qf,1)/R; %probability of a spike on a channel
%     q=1-p; %probability of no spike on a channel
%     pexp=zeros(npat,1);
%     spike=zeros(1,C); nospike=zeros(1,C);
%     %for j=1:length(desort); jj=ind(j); %time index of each unique pattern
%     for j=1:npat
%         jj=ind1(j); %time index of each unique pattern
%         for i=1:C; spike(i)=p(i).^qf(jj,i); end %spikes on channels
%         for i=1:C; nospike(i)=q(i).^(1-qf(jj,i)); end %no spike on channels
%         pexp(j)=prod(spike)*prod(nospike); %expected prob of observed pattern
%         %assuming independence of channels
%     end
%     nexp=pexp*R; %number of expected events
%     pspu=psp(count1>nexp); %prob of unitary events
%     npspu=length(pspu); %number of unitary events
%     indu=ind1(count1>nexp); %index of unitary events
%     qu=qf(indu,:); %unitary spike patterns

 %------entropy for subsets used for complexity-----------------------------
    sumE=0;
    for i=1:C
        clear M;
        qfr=qfs; qfr(:,i)=[]; %delete one channel
        M=cellstr(qfr); %cell array of spike patterns as strings
        [~, ind2]=unique(sort(M)); 
        count2=diff([0;ind2]); %counts of unique patterns
        count2(1)=[]; %remove no spike pattern ********************************
        pspr=count2/sum(count2); %probability of unique spike patterns
        entros=-sum(pspr.*log2(pspr)); %entropy of subset reduced by one channel
        sumE=sumE+entros; %sum up subset entropies
        clear ind2 pspr;
    end   
    
%-------surrogate data-------------------------------------------
        qs=qf; %qs(sum(qs,2)==0,:)=[]; %remove no-spike patterns - use for EI only
        nr=1; %randomization repeats     
        
        %HAVE TO CHANGE THIS BACK FOR GOOD RESULTS
        
        sumE2_temp=zeros(nr,1);
        for i=1:nr %repeat randomization nr times, 10 times is good
        
%             %shuffle spike trains
%             maxshift=300; %make this pre/post duration*2
%             shift=randperm(maxshift)-ceil(maxshift/2);
%             shift(shift==0)=[];
%             for j=1:C     
%                 qs(:,j)=circshift(qs(:,j),shift(j)); 
%             end
        
        %random spike generation
        spr=sum(qs)/length(qs); %spike rate matrix
        for j=1:C
            qs(:,j)=binornd(1,spr(j),length(qs),1);
            %for k=1:R;qs(k,j)=rand<spr(j);end %alternative calculation
        end
         
        %figure(3);imagesc(qs',[-1 1]);
        qfs=num2str(qs,2); M=cellstr(qfs); %cell array of spike patterns as strings
        [~, ind3]=unique(sort(M));
        counts=diff([0;ind3]); %counts of unique patterns
        psps=counts/sum(counts); %probability of unique spike patterns, surrogate  
        entro2_tmp(i)=-sum(psps.*log2(psps)); %from surrogate data
        %clear M qfs ind3 counts psps;
    
        %subsets for shuffled complexity CI    
        for j=1:C
            clear M;
            qfr=qfs; qfr(:,j)=[]; %delete one channel
            M=cellstr(qfr); %cell array of spike patterns as strings
            [~, ind2]=unique(sort(M)); 
            count2=diff([0;ind2]); %counts of unique patterns            
            count2(1)=[]; %remove no spike pattern ********************************           
            pspr=count2/sum(count2); %probability of unique spike patterns
            entros=-sum(pspr.*log2(pspr)); %entropy of subset reduced by one channel
            sumE2_temp(i)=sumE2_temp(i)+entros; %sum up subset entropies
            clear ind2 pspr;
        end
        sumE2=mean(sumE2_temp);
        
    end
    npats=length(counts); %number of surrogate spike patterns (from the last run)
        
    %entropies:
%     entro0=-sum(pexp.*log2(pexp)); %from expected probabilities
    entro1=-sum(psp.*log2(psp));   %from unique spike patterns #probably HX
    entro2=mean(entro2_tmp); %from surrogate data 
    
    %hey="hey"
    %size(psp)
    try
    psp(1)=[];
    entro3=-sum(psp.*log2(psp));   %from unique spike patterns
    
    p=sum(qf,1)/R; q=1-p; %probability of a spike or no spike on a channel
    entro4=-sum(p.*log2(p))-sum(q.*log2(q)); %used for total correlation # this is my H_x_i(i)
    
    TC=entro4-entro1; %total correlation; independent minus unique entropy # so this is my integration
    EI=entro2-entro3; %interaction entropy; surrogate minus unique NO ZEROS
    CI=sumE-(C-1)*entro1; %complexity 
    TC2=entro4-entro2; %from shuffled data
    CI2=sumE2-(C-1)*entro2; %from shuffled data
catch
        TC=0; EI=0; CI=0; psp=0;
        ALLEn=[0, 0, 0, 0, 0,0,0,0];
    end
    ALLEn=[npat,npats,TC,TC2,EI,CI,CI2,entro1];
end
