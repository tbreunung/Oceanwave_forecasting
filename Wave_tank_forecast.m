clear all
close all
clc

src='HPPW_f0_1p30_ak_0p14_sbFrac_0p05_unb_0p00_psi_1p00_Omegam_1p00_eta.mat';
%src='J_RAO_TP_0p80_Hs_14_foc_10_nbT_1p60_gamma_3p30_eta.mat';
load(src)

%%

station=1;

Fs=400;
T=1/Fs;
shift=1000;
N_pts=1000;
N_wdw=floor((length(wh(station,:))-N_pts)/shift);
%ampls_single_sided=zeros(N_wdw,N_pts/2+1);
figure
% freqs=Fs*(0:(N_pts/2))/N_pts;
% for iter_wdw=1:N_wdw
% ampls_two_sided=abs(fft(wh(station,(iter_wdw-1)*shift+(1:N_pts))./N_pts));
% ampls_single_sided=ampls_two_sided(1:N_pts/2+1);
% ampls_single_sided(2:end)=2*ampls_single_sided(2:end);
% 
% plot(freqs,ampls_single_sided)
% hold on
% end

N_pts=length(wh(station,:));
 
freqs=Fs*(0:(N_pts/2))/N_pts;
ampls_two_sided=abs(fft(wh(station,:))./N_pts);
ampls_single_sided=ampls_two_sided(1:N_pts/2+1);
ampls_single_sided(2:end)=2*ampls_single_sided(2:end);

plot(freqs,ampls_single_sided)
xlabel('Frequency in Hz')
ylabel('Energy density')


 
xlim([0 5])

%%
eps=1*10^-3;
alpha=10;

N_pts_per_wave=3;
N_windows=1;
N_MC=1000;
N_fmodes=1;

err_const=NaN(N_windows,1);
err0=NaN(N_windows,1);
err_MC=NaN(N_windows,1);
cost_const=NaN(N_windows,1);
cost_SF=NaN(N_windows,1);
N_waves=[60 60];

samplerate=400;


%err_succ=NaN(N_windows,window_size(2));

 for jj=1:N_windows
   
   
     w_start=  1;%randi(info.Size-sum(window_size));%41196265;%6781226;%10398650; %9720175;%
     
     hz=wh(station,:);%ncread(src,'xyzZDisplacement',w_start,window_size).';
     
     [tmp_wave_height ,tmp_wave_idx,zero_idx, ~]=my_wave_height_filter(hz ,N_pts_per_wave);
     
      
     
     ts=[tmp_wave_idx ];%  tmp_t0];%[1:window_size];%
     hs=[tmp_wave_height];%   zeros(1,length(tmp_t0))];%hz(jj,:);%
     
     [ts,id]=sort(ts);
     hs=hs(id);
    
      
     t1=ts(1:N_waves(1));
     h1=hs(1:N_waves(1));
    
     if N_fmodes<9
         str_fmode=append('sin',num2str(N_fmodes));
         my_fit=fit(t1.',(h1-mean(h1)).',str_fmode);
         err_const(jj)=sum((h1-mean(h1)-my_fit(t1).').^2);
         
         coeff2=coeffvalues(my_fit);
         ws=sort(coeff2(2:3:end));
         tic
         [ws_MC, err_MC(jj)]= my_Freq_MC_fit(t1,h1,err_const(jj),N_fmodes,ws,N_MC, 1/Fs);
         toc
     else
         ws=zeros(1,N_fmodes);
         err_const(jj)=sum((h1-mean(h1)).^2);
         tic
         [ws_MC, err_MC(jj)]= my_Freq_MC_fit2(t1,h1,err_const(jj),N_fmodes,ws,N_MC,1/Fs);
         toc
     end
     N_waves(2)=length(hs)-1-N_waves(1);
     hz_slowflow=cell(N_windows,N_waves(2)+1);
    slow_vars=cell(N_windows,N_waves(2)+1);%NaN(N_windows,N_waves(2)+1,2*N_fmodes+1,N_waves(1)+N_waves(2));
    err_SF=NaN(N_windows,N_waves(2)+1);
     [hz_slowflow{jj,1} ,slow_vars{jj,1}, err_SF(jj,1)]= my_SF_fit(t1, h1,ws_MC,eps,alpha);
%      =svars_tmp;
        
%     =err_tmp;
     t_tmp=cell(N_waves(2)+1,1);
     h_tmp=cell(N_waves(2)+1,1);%,N_waves(1)+N_waves(2));
     h_tmp{1}=h1;
     t_tmp{1}=t1;
     tic
     parfor ii=1:length(ts)-N_waves(1)-1
         t_tmp{ii+1}= ts(1:N_waves(1)+ii);
         h_tmp{ii+1}= hs(1:N_waves(1)+ii);%[h_tmp(ii,2:end) hs(N_waves(1)+ii)];
         [hz_slowflow{jj,ii+1} ,slow_vars{jj,ii+1},err_SF(jj,ii+1)]= my_SF_fit(t_tmp{ii+1}, h_tmp{ii+1},ws_MC,eps,alpha);
       

     end
     toc
       [hz_slowflow_succ, slow_vars_succ,err_succ]= my_SF_fit_succ(ts(N_waves(1)+1:end), hs(N_waves(1)+1:end),ws_MC,slow_vars{1,1}(:,end),t_tmp{1}(end),eps);
%      hz2=ncread(src,'xyzZDisplacement',w_start+zero_idx(end)-1,window_size(2)).';
%      
%      [tmp_wave_height ,tmp_wave_idx,~, ~]=my_wave_height_filter(hz2 ,N_pts_per_wave);
%      
%      
%      
%      t2=[tmp_wave_idx ];%  tmp_t0];%[1:window_size];%
%      h2=[tmp_wave_height];%   zeros(1,length(tmp_t0))];%hz(jj,:);%
%      
%      [t2,id]=sort(t2);
%      h2=h2(id);
%      
%      [hz_slowflow2,slow_vars2,err_SF2]= my_SF_fit_succ(zero_idx(end)+t2-1, h2,ws_MC,slow_vars1(:,end),t1(end),eps);
%      err_SF(jj)=err_SF1+err_SF2;
%      
%      err0(jj)=sum(h1.^2)+sum(h2.^2);
%      tmp=NaN(1,window_size(2));
%      for ll=1:length(t2)
%          tmp(ll)=(norm(h2(1:ll)-hz_slowflow2(1:ll)))/norm(h2(1:ll));
%      end
%      err_succ(jj,:)=tmp;
    
      jj
 end
 
%  hz3=ncread(src,'xyzZDisplacement',w_start,window_size(1)+window_size(2)).';
%  
%  [h3 ,t3,~, ~]=my_wave_height_filter(hz3 ,N_pts_per_wave);
  

%%


err_0=cell2mat(cellfun(@norm,h_tmp,'Un',0));
rel_err=sqrt(err_SF).'./err_0;

t_end=cell2mat(cellfun(@(a) a(end),t_tmp,'Un',0));
figure
plot(t_end./samplerate,rel_err)
xlabel('time in s')
ylabel('realtive least squares error')


window_idx=1;

figure
%subplot(2,1,1)
plot((1:length(hz))./samplerate,hz)
hold on
hz_fit=SF2fulltime(t_tmp{window_idx}(1):t_tmp{window_idx}(end),t_tmp{window_idx},slow_vars{window_idx},ws_MC,N_fmodes);
%hz_fit_succ=SF2fulltime(ts(N_waves(1)+1):ts(end),ts(N_waves(1)+1:end),slow_vars_succ,ws_MC,N_fmodes);


plot((t_tmp{window_idx}(1):t_tmp{window_idx}(end))./samplerate,hz_fit)

%plot(ts(N_waves(1)+1):ts(end),hz_fit_succ,'-.')

plot(t_tmp{window_idx}./samplerate,h_tmp{window_idx},'xr')

plot(t_tmp{window_idx}./samplerate,hz_slowflow{1,window_idx},'s')
legend('data','fit')
xlim([t_tmp{window_idx}(1)./samplerate t_tmp{window_idx}(end)./samplerate])
ylabel('surface elevation (mm)')
xlabel('time (s)')
% subplot(2,1,2)
% plot(slow_vars{window_idx}.')

%%
N_dim=2*N_fmodes+1;

 %fit_steps=99;
%N_waves(2)-fit_steps;
delay_dim=50;

fit_window=1;
lag_steps= 1;
mode='DMD';
sample_data=10;
slow_vars_fit=slow_vars{fit_window};
forecast_steps=800;%
 

tt=t_tmp{fit_window}(1):sample_data:t_tmp{fit_window}(end)-lag_steps;
% for jj=1:N_dim
slow_vars_interp =interp1(t_tmp{fit_window},slow_vars_fit.',tt,'spline','extrap').';
% end
 

m0=mean(slow_vars_interp,2);
m1=std(slow_vars_interp,0,2)./std(slow_vars_interp,0,2);
slow_vars_nomean=(slow_vars_interp-m0)./m1;



%[lpc_coeffs,~]=lpc(slow_vars_nomean.',delay_dim);
% ar_coeffs=zeros(N_dim,delay_dim+1);
% es=zeros(N_dim,1);
% refl_coeffs=zeros(N_dim,delay_dim);
% for iter_dim=1:N_dim
% [r,lg] = xcorr(slow_vars_nomean(iter_dim,:).','biased');
% r(lg<0) = [];
% 
% [ar_coeffs(iter_dim,:), es(iter_dim) ,refl_coeffs(iter_dim,:)]=levinson(r,delay_dim);
% end
 A_mat_DMD=zeros(N_dim,delay_dim*N_dim);
 for jj=1:N_dim
     zs_mat=[];
     for dd=1:delay_dim
         zs_mat=[zs_mat; slow_vars_nomean(jj,dd:end-delay_dim+dd)];
     end
     V1=zs_mat(:,1:end-1);
     V2=zs_mat(:,2:end);
     [U,S,V]=svd(V1,'econ');
     svds=1./diag(S);
     svds(svds>1)=0;
     A_tmp=V2*V*diag(svds)*U';
     A_mat_DMD(jj,jj:N_dim:end)=A_tmp(end,:);
 end
 
 A_mat_L2=zeros(N_dim,delay_dim*N_dim);
 for jj=1:N_dim
     zs_mat=[];
     for dd=1:delay_dim
         zs_mat=[zs_mat; slow_vars_nomean(jj,dd:end-delay_dim+dd-1)];
     end
     %reshape(cell2mat(cellfun(@(a) a(jj,:),zs_nomean,'UniformOutput',false)),delay_dim,[]).';%squeeze(zs(:,jj,:));
     %zs_mat(fit_steps+1:end,:)=[];
     
     ys_vec=slow_vars_nomean(jj,delay_dim+1:end).';%-slow_vars_nomean(jj,delay_dim:end-1).';%ys_delta(jj,:).';%
     A_vec=zs_mat.'\ys_vec;
     A_mat_L2(jj,jj:N_dim:end)=A_vec;
     
 end
  

    
  
 
 
 
 

%%

numResponses= N_dim;%
numFeatures= N_dim;%
numHiddenUnits=32;


layers = [ ...
    sequenceInputLayer(numFeatures)
    lstmLayer(numHiddenUnits,'OutputMode','sequence')
     fullyConnectedLayer(200)
    dropoutLayer(0.1)
     lstmLayer(numHiddenUnits/2,'OutputMode','sequence')
     fullyConnectedLayer(80)
     dropoutLayer(0.1)
     lstmLayer(numHiddenUnits/4,'OutputMode','sequence')
     
    % lstmLayer(numHiddenUnits/8,'OutputMode','sequence')
     
   fullyConnectedLayer(80)
    dropoutLayer(0.1)
%      lstmLayer(numHiddenUnits,'OutputMode','sequence')
%     fullyConnectedLayer(50)
%     dropoutLayer(0.1)
     fullyConnectedLayer(numResponses)
    regressionLayer];


options = trainingOptions('adam', ...
    'MaxEpochs',150, ...
    'GradientThreshold',1, ...
    'InitialLearnRate',0.005, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropPeriod',125, ...
    'LearnRateDropFactor',0.2, ...
    'Verbose',0, ...
    'Plots','training-progress');


% nets=cell(N_dim,1);
%  parfor iter_dim=1:N_dim
%  nets{iter_dim}    = trainNetwork(slow_vars_nomean(iter_dim,1:end-1) ,slow_vars_nomean(iter_dim,2:end) ,layers,options);
%  end
 net  = trainNetwork(slow_vars_nomean(:,1:end-1) ,slow_vars_nomean(:,2:end) ,layers,options);

%%
N_t=length(slow_vars_nomean(1,:));
slow_vars_forecast_L2=zeros(N_dim,N_t+forecast_steps);
slow_vars_forecast_DMD=zeros(N_dim,N_t+forecast_steps);
slow_vars_forecast_triv=zeros(N_dim,N_t+forecast_steps);
slow_vars_forecast_LSTM=zeros(N_dim,N_t+forecast_steps);

slow_vars_forecast_L2(:,1:N_t)=slow_vars_nomean;
slow_vars_forecast_DMD(:,1:N_t)=slow_vars_nomean;
slow_vars_forecast_triv(:,1:N_t)=slow_vars_nomean;
slow_vars_forecast_LSTM(:,1:N_t)=slow_vars_nomean;



 for jj=1:forecast_steps
    
    tmp=slow_vars_forecast_L2(:,N_t+jj+(-delay_dim:-1));
    slow_vars_forecast_L2(:,N_t+jj)=A_mat_L2*tmp(:);%slow_vars_forecast_own(:,jj+delay_dim-1)+
    
    tmp=slow_vars_forecast_DMD(:,N_t+jj+(-delay_dim:-1));
    slow_vars_forecast_DMD(:,N_t+jj)=A_mat_DMD*tmp(:);%slow_vars_forecast_own(:,jj+delay_dim-1)+
    
     
    %   pred=zeros(N_dim,N_t+jj-1);
    %for ii=1:N_dim
      
     %[nets{ii},pred(ii,:)]=predictAndUpdateState(nets{ii},slow_vars_forecast_LSTM(ii,1:N_t+jj-1),'MiniBatchSize',1);
   % end
    [ net, pred] =predictAndUpdateState(net ,slow_vars_forecast_LSTM(:,1:N_t+jj-1),'MiniBatchSize',1);   
    slow_vars_forecast_LSTM(:,N_t+jj)= pred(:,end) ; 

    slow_vars_forecast_triv(:,N_t+jj)=slow_vars_nomean(:,end);

end
slow_vars_forecast_L2=m0+slow_vars_forecast_L2.*m1;
slow_vars_forecast_DMD=m0+slow_vars_forecast_DMD.*m1;
slow_vars_forecast_LSTM=m0+slow_vars_forecast_LSTM.*m1;
slow_vars_forecast_triv=m0+slow_vars_forecast_triv.*m1;

slow_vars_forecast_DMD(:,1:N_t)=[];
slow_vars_forecast_L2(:,1:N_t)=[];
slow_vars_forecast_LSTM(:,1:N_t)=[];
slow_vars_forecast_triv(:,1:N_t)=[];
slow_vars_forecast_mean=repmat(m0,1,forecast_steps);
%slow_vars_forecast_lpc(:,1:delay_dim)=[];

%%
figure
%id=randi(N_dim);
window_idx=find(cell2mat(cellfun(@(a)a(end),t_tmp,'Un',0))>t_tmp{fit_window}(end)+forecast_steps*sample_data-lag_steps,1,'first');

tmp=plot(t_tmp{fit_window}./samplerate,slow_vars{fit_window} );
pl(1,:)=tmp(1,:);
%plot(slow_vars2{fit_window}.')
hold on
%set(gca,'ColorOrderIndex',1)
%tmp=plot(tt./samplerate,slow_vars_interp ,'-.');
%pl(1,:)=tmp(1,:);
set(gca,'ColorOrderIndex',1)
tmp=plot(t_tmp{window_idx}./samplerate,slow_vars{window_idx} ,'--');
pl(2,:)=tmp(1,:);
%plot(slow_vars2{window_idx}.','--')
% set(gca,'ColorOrderIndex',1)
% plot(((1-delay_dim:forecast_steps)+t_tmp{fit_window}(end)-lag_steps)./samplerate,slow_vars_forecast_lpc,'x-')





set(gca,'ColorOrderIndex',1)
tmp=plot(((1:forecast_steps)*sample_data+t_tmp{fit_window}(end)-lag_steps)./samplerate,slow_vars_forecast_L2 ,'-d','MarkerIndices',1:10:length(slow_vars_forecast_L2(1,:)));
pl(3,:)=tmp(1,:);
set(gca,'ColorOrderIndex',1)
tmp=plot(((1:forecast_steps)*sample_data+t_tmp{fit_window}(end)-lag_steps)./samplerate,slow_vars_forecast_DMD ,'-s','MarkerIndices',1:10:length(slow_vars_forecast_L2(1,:)));
pl(4,:)=tmp(1,:);
set(gca,'ColorOrderIndex',1)
tmp=plot(((1:forecast_steps)*sample_data+t_tmp{fit_window}(end)-lag_steps)./samplerate,slow_vars_forecast_LSTM,'-x','MarkerIndices',1:10:length(slow_vars_forecast_L2(1,:)));
pl(5,:)=tmp(1,:);

plot(t_tmp{fit_window}(end)/samplerate.*[1 1],ylim,'--k')


xlabel('time in s')
ylabel('amplitudes a_j')
legend(pl,'Measurements for fit','true time series','L2-forecast','DMD-forecast','LSTMs','location','NorthWest' )

%ylim([-1 1])

slow_vars_truth=zeros(1,forecast_steps);
err_DMD_SF=zeros(1,forecast_steps);
err_L2_SF=zeros(1,forecast_steps);
err_triv_SF=zeros(1,forecast_steps);
%err_arma_SF=zeros(1,forecast_steps);
err_LSTM_SF=zeros(1,forecast_steps);
err_mean_SF=zeros(1,forecast_steps);



for iter_forecast=1:forecast_steps
    %     if iter_forecast<lag_steps
    %         tt=t_tmp{fit_window}(end)-lag_steps+iter_forecast;
    %         for jj=1:N_dim
    %             slow_vars_truth(jj,iter_forecast)=interp1(t_tmp{fit_window},slow_vars_fit(jj,:),tt,'spline','extrap');
    %         end
    %     else
    if iter_forecast<lag_steps+1
        for jj=1:N_dim
            tmp=slow_vars{fit_window};
            slow_vars_truth(jj,1:iter_forecast)=interp1(t_tmp{fit_window},tmp(jj,:),(t_tmp{fit_window}(end)-lag_steps)+(1:iter_forecast),'spline','extrap');
            
            %slow_vars_truth(jj,1:iter_forecast)= tmp(jj,(t_tmp{fit_window}(end)-lag_steps)+(1:iter_forecast));
            %slow_vars_forecast_triv(jj,delay_dim+(1:iter_forecast))=slow_vars_truth(jj,1:iter_forecast);
        end
        
    else
        t_idx=t_tmp{fit_window}(end)+iter_forecast*sample_data-lag_steps;
        window_idx=find(cell2mat(cellfun(@(a)a(end),t_tmp,'Un',0))>t_idx,1,'first');
        for jj=1:N_dim
            tmp=slow_vars{window_idx};
            slow_vars_truth(jj,1:iter_forecast)=interp1(t_tmp{window_idx},tmp(jj,:),(t_tmp{fit_window}(end)-lag_steps+1:sample_data:t_idx),'spline','extrap');

            %slow_vars_truth(jj,1:iter_forecast)=tmp(jj,(t_tmp{fit_window}(end)-lag_steps)+(1:iter_forecast));

        end
    end
    %end
    err_DMD_SF(iter_forecast)=mean(vecnorm(slow_vars_forecast_DMD(:,1:iter_forecast)-slow_vars_truth(:,1:iter_forecast)).^2);
    err_L2_SF(iter_forecast)=mean(vecnorm(slow_vars_forecast_L2(:,1:iter_forecast)-slow_vars_truth(:,1:iter_forecast)).^2);
    err_triv_SF(iter_forecast)=mean(vecnorm(slow_vars_forecast_triv(:,1:iter_forecast)-slow_vars_truth(:,1:iter_forecast)).^2);
    % err_arma_SF(iter_forecast)=mean(vecnorm(slow_vars_forecast_ARMA(:,1:iter_forecast)-slow_vars_truth(:,1:iter_forecast)).^2);
    err_LSTM_SF(iter_forecast)=mean(vecnorm(slow_vars_forecast_LSTM(:,1:iter_forecast)-slow_vars_truth(:,1:iter_forecast)).^2);
    err_mean_SF(iter_forecast)=mean(vecnorm(slow_vars_forecast_mean(:,1:iter_forecast)-slow_vars_truth(:,1:iter_forecast)).^2);
    
    
end
%%
figure
% for iter_idx=1:min([N_dim 12])
% subplot(6,4,iter_idx)
% pl_tmp=plot((1:forecast_steps)*sample_data./samplerate,slow_vars_forecast_lpc(iter_idx,:));
% pl(1,:)=pl_tmp(1);
% hold on
% %set(gca,'ColorOrderIndex',1)
% pl_tmp=plot((1:forecast_steps)*sample_data./samplerate,slow_vars_forecast_own(iter_idx,:));
% pl(2,:)=pl_tmp(1,:);
% 
% pl_tmp=plot((1:forecast_steps)*sample_data./samplerate,slow_vars_forecast_LSTM(iter_idx,:));
% 
% %set(gca,'ColorOrderIndex',1)
% pl_tmp=plot((1:forecast_steps)*sample_data./samplerate,slow_vars_truth(iter_idx,:));
% pl(3,:)=pl_tmp(1,:);
% 
% %plot((1:forecast_steps)./samplerate,slow_vars_forecast_ARMA(iter_idx,:));
% xlabel('time in s')
% ylabel(['amplitudes a_' num2str(iter_idx)])
% %legend(pl,'LPC','L2-Fit','truth')
% %ylim([-1 1])
% 
% end
% subplot(6,4, 13:24 )
plot(((1:forecast_steps)*sample_data+t_tmp{fit_window}(end)-lag_steps)./samplerate,err_L2_SF)
hold on
plot(((1:forecast_steps)*sample_data+t_tmp{fit_window}(end)-lag_steps)./samplerate,err_DMD_SF)
plot(((1:forecast_steps)*sample_data+t_tmp{fit_window}(end)-lag_steps)./samplerate,err_LSTM_SF)
plot(((1:forecast_steps)*sample_data+t_tmp{fit_window}(end)-lag_steps)./samplerate,err_triv_SF)
plot(((1:forecast_steps)*sample_data+t_tmp{fit_window}(end)-lag_steps)./samplerate,err_mean_SF)
 %plot((1:forecast_steps)./samplerate,err_arma_SF)

xlabel('time in s')
ylabel('mean square error slow flow')
lg=legend('L2-fit','DMD','LSTM','last sample','training mean');
set(lg,'Location','Northwest')

ylim([0 max([err_triv_SF])])

xlim([(t_tmp{fit_window}(end)-lag_steps)./samplerate (forecast_steps*sample_data+t_tmp{fit_window}(end)-lag_steps)./samplerate])
figure

window_idx=find(cell2mat(cellfun(@(a)a(end),t_tmp,'Un',0))>t_tmp{fit_window}(end)+forecast_steps*sample_data-lag_steps,1,'first');
hz_slowflow_full_truth=SF2fulltime(t_tmp{fit_window}(end)-lag_steps+(1:forecast_steps*sample_data),t_tmp{window_idx},slow_vars{window_idx},ws_MC,N_fmodes);
hz_slowflow_full_forecast=SF2fulltime(t_tmp{fit_window}(end)-lag_steps+(1:forecast_steps*sample_data),t_tmp{fit_window}(end)-lag_steps+(1:forecast_steps)*sample_data,slow_vars_forecast_LSTM,ws_MC,N_fmodes);

plot((1:length(hz))./samplerate,hz)
hold on

plot((t_tmp{fit_window}(end)-lag_steps+(1:forecast_steps*sample_data))./samplerate,hz_slowflow_full_truth)
plot((t_tmp{fit_window}(end)-lag_steps+(1:forecast_steps*sample_data))./samplerate,hz_slowflow_full_forecast)
xlim((t_tmp{fit_window}(end)-lag_steps+[1 forecast_steps*sample_data])./samplerate)
 
idx_start=find(t_tmp{window_idx}>(t_tmp{fit_window}(end)),1,'first');%

set(gca,'ColorOrderIndex',1)
plot((t_tmp{window_idx})./samplerate,h_tmp{window_idx},'x')
%plot((t_tmp{window_idx}(idx_start:end-1))./samplerate,hz_slowflow_full_forecast(t_tmp{window_idx}(idx_start:end-1)-t_tmp{fit_window}(end)+lag_steps),'dg')


[fit_wave_height ,fit_time,~, ~]=my_wave_height_filter(hz_slowflow_full_truth ,N_pts_per_wave);
plot((fit_time+t_tmp{fit_window}(end)-lag_steps)./samplerate,fit_wave_height,'s')
[forecast_wave_height ,forecast_time,~, ~]=my_wave_height_filter(hz_slowflow_full_forecast ,N_pts_per_wave);
plot((forecast_time+t_tmp{fit_window}(end)-lag_steps)./samplerate,forecast_wave_height,'d')



legend('full data','fit','forecast')
xlabel('time  (s)')
ylabel('surface elevation (mm)')

forecast_time=forecast_time+t_tmp{fit_window}(end)-lag_steps;
fit_time=fit_time+t_tmp{fit_window}(end)-lag_steps;
for iter_waves=idx_start:length(h_tmp{window_idx})-1
    
    [~,idx]=min(abs(t_tmp{window_idx}(iter_waves)-forecast_time));
    if h_tmp{window_idx}(iter_waves)*forecast_wave_height(idx)>0
       h_forecast(iter_waves+1-idx_start)= forecast_wave_height(idx);
    else
        if forecast_time(idx)<t_tmp{window_idx}(iter_waves)
            if   idx<length(forecast_time)
                h_forecast(iter_waves+1-idx_start)= forecast_wave_height(idx+1);
            else
                h_forecast(iter_waves+1-idx_start)= forecast_wave_height(idx-1);
            end
        else
             if   idx>1
                h_forecast(iter_waves+1-idx_start)= forecast_wave_height(idx-1);
             else
                h_forecast(iter_waves+1-idx_start)= forecast_wave_height(idx+1);
             end
   
        end
    end
     [~,idx]=min(abs(t_tmp{window_idx}(iter_waves)-fit_time));
    if h_tmp{window_idx}(iter_waves)*fit_wave_height(idx)>0
       h_fit(iter_waves+1-idx_start)= fit_wave_height(idx);
    else
        if fit_time(idx)<t_tmp{window_idx}(iter_waves) 
            if idx<length(fit_time)
                h_fit(iter_waves+1-idx_start)= fit_wave_height(idx+1);
            else
                h_fit(iter_waves+1-idx_start)= fit_wave_height(idx-1);    
            end
        else
            if idx>1
            h_fit(iter_waves+1-idx_start)= fit_wave_height(idx-1);
            else
                h_fit(iter_waves+1-idx_start)= fit_wave_height(idx+1);
            end
            
        end
    end
    
end



err_forecast=mean(abs(h_tmp{window_idx}(idx_start:end-1)-h_forecast))

err_fit=mean(abs(h_tmp{window_idx}(idx_start:end-1)-h_fit))

err_gauss=mean(abs(abs(h_tmp{window_idx}(idx_start:end-1))-mean(abs(h_tmp{fit_window}))))
%corr(forecast_wave_height.',h_tmp{window_idx}(idx_start:end-1).')

figure
plot(abs(h_tmp{window_idx}(idx_start:end-1)),'x')
hold on
plot(abs(h_fit),'o')
plot(abs(h_forecast),'s')

 
%%
function [wave_height,wave_idx,zero_crossing,t0]=my_wave_height_filter(hz,N_pts_per_wave)

hz_mean=mean(hz);
zero_crossing=find(abs(diff(sign(hz-hz_mean) ))==2);

insig_waves =diff(zero_crossing)<N_pts_per_wave;

zero_crossing( insig_waves )=[];



wave_height=zeros(1,length(zero_crossing)-1);
wave_idx=zeros(1,length(zero_crossing)-1);

for ii=1:length(zero_crossing)-1
        [~, tmp_idx]=max(abs(hz((zero_crossing(ii)+1):(zero_crossing(ii+1)))));
        wave_idx(ii)=tmp_idx+zero_crossing(ii);
        wave_height(ii)= hz(wave_idx(ii));
        
end


idx=find(diff(sign(wave_height))==0,1);
while isempty(idx)~=1
    tmp_heights=[wave_height(idx) wave_height(idx+1)];
    tmp_wave_idxs=[wave_idx(idx) wave_idx(idx+1)];
    [~, tmp_idx]=max(abs( tmp_heights));
    wave_height(idx)=tmp_heights( tmp_idx);
    wave_idx(idx)=tmp_wave_idxs(tmp_idx);
    wave_height(idx+1)=[];
    wave_idx(idx+1)=[];
    zero_crossing(idx+1)=[];
    idx=find(diff(sign(wave_height))==0,1);
end

for ii=1:length(zero_crossing)
   % if zero_crossing(ii+1)-zero_crossing(ii)>N_pts_per_wave-1
        t0(ii)=zero_crossing(ii)-hz(zero_crossing(ii))/(hz(zero_crossing(ii)+1)-hz(zero_crossing(ii)));
end
end

function [hz_slowflow, slow_vars,err]= my_SF_fit(ts, hs,ws,eps,alpha)
N_t=length(ts);
N_fmodes=length(ws);

M1=zeros(N_t,(2*N_fmodes+1)*N_t);
M2=zeros((2*N_fmodes+1)*N_t,(2*N_fmodes+1)*N_t);
M3=zeros(N_t,(2*N_fmodes+1)*N_t);

for tt=1:N_t
    idx0=1+(tt-1)*(2*N_fmodes+1);
    idx1=idx0+N_fmodes-1;
    idx2=idx1+N_fmodes;
    M1(tt,idx0)=1;
    M1(tt,idx0+1:idx1+1)=  cos(ws.*ts(tt));
    M1(tt,idx1+2:idx2+1)= sin(ws.*ts(tt));
    M3(tt,idx0+1:idx1+1)= -ws.*sin(ws.*ts(tt));
    M3(tt,idx1+2:idx2+1)= ws.*cos(ws.*ts(tt));
    
    if tt<N_t
        del_t=(ts(tt+1)-ts(tt));
        idx1=1+(tt-1)*(2*N_fmodes+1);
        M2(idx1,idx1)=-1/del_t;
        M2(idx1,idx1+2*N_fmodes+1)=1/del_t;
        %M2(idx1,idx1-2*N_fmodes-1)=1/del_t;
        for iter_fmodes=1:N_fmodes
            M2(idx1+iter_fmodes,idx1+iter_fmodes)=-1/del_t;
            M2(idx1+iter_fmodes,idx1+iter_fmodes+2*N_fmodes+1)=1/del_t;
            M2(idx1+iter_fmodes+N_fmodes,idx1+iter_fmodes+N_fmodes)=-1/del_t;
            M2(idx1+iter_fmodes+N_fmodes,idx1+iter_fmodes+3*N_fmodes+1)=1/del_t;
        end
    else
%         del_t=(ts(tt)-ts(tt-1));
%         idx1=1+(N_t-1)*(2*N_fmodes+1);
% 
%         M2(idx1,idx1)=1/del_t;
%         M2(idx1,idx1-(2*N_fmodes+1))=M2(idx1,idx1-(2*N_fmodes+1))-1/del_t;
% 
%         for iter_fmodes=1:N_fmodes
%            
%             
%             M2(idx1+iter_fmodes,idx1+iter_fmodes)= 1/del_t;
%             M2(idx1+iter_fmodes,idx1+iter_fmodes-(2*N_fmodes+1))= M2(idx1+iter_fmodes,idx1+iter_fmodes-(2*N_fmodes+1))-1/del_t;
%             
%             M2(idx1+iter_fmodes+N_fmodes,idx1+iter_fmodes+N_fmodes )=1/del_t;
%             M2(idx1+iter_fmodes+N_fmodes,idx1+iter_fmodes-N_fmodes-1)=M2(idx1+iter_fmodes+N_fmodes,idx1+iter_fmodes-N_fmodes-1)-1/del_t;
%             
%  
%         end
    end
end
 
Mbig=(M1.'*M1+1/eps.*M2.'*M2+alpha.*M3.'*M3);
Vs=Mbig\(M1.'*hs.');

 
hz_slowflow=zeros(1,N_t);
slow_vars=zeros(2*N_fmodes,N_t);
for tt=1:N_t
    idx1=1+(tt-1)*(2*N_fmodes+1);
    
     hz_slowflow(tt)=Vs(idx1);
    slow_vars(1,tt)=Vs(idx1);
    for iter_fmodes=1:N_fmodes
        hz_slowflow(tt)=hz_slowflow(tt)+Vs(idx1+iter_fmodes)*cos(ws(iter_fmodes)*ts(tt))...
            + Vs(idx1+iter_fmodes+N_fmodes)*sin(ws(iter_fmodes)*ts(tt));
    slow_vars(1+iter_fmodes,tt)=Vs(idx1+iter_fmodes);        
    slow_vars(1+N_fmodes+iter_fmodes,tt)=Vs(idx1+N_fmodes+iter_fmodes);    
    end
end
 

err=sum((hz_slowflow-hs).^2);
 

end


function [hz_slowflow, slow_vars,err]= my_SF_fit2(ts, hs,ws,eps)
N_t=length(ts);
N_fmodes=length(ws);

M1=zeros(N_t,(2*N_fmodes)*N_t);
M2=zeros((2*N_fmodes)*N_t,(2*N_fmodes)*N_t);
M3=zeros(N_t,(2*N_fmodes)*N_t);

for tt=1:N_t
    idx0=1+(tt-1)*(2*N_fmodes);
    idx1=idx0+N_fmodes-1;
    idx2=idx1+N_fmodes;
    M1(tt,idx0:idx1)= abs( cos(ws.*ts(tt)));
    M1(tt,(idx1+1):idx2)= abs(sin(ws.*ts(tt)));
    if cos(ws.*ts(tt))>0
        M3(tt,idx0:idx1)= -ws.*sin(ws.*ts(tt));
    else
        M3(tt,idx0:idx1)= ws.*sin(ws.*ts(tt));
    end 
    if sin(ws.*ts(tt))>0
        M3(tt,idx1+1:idx2)= ws.*cos(ws.*ts(tt));
    else
         M3(tt,idx1+1:idx2)= -ws.*cos(ws.*ts(tt));
    end
    
    if tt<N_t
        del_t=(ts(tt+1)-ts(tt));
        idx1=(tt-1)*(2*N_fmodes);
        %M2(idx1,idx1-2*N_fmodes-1)=1/del_t;
        for iter_fmodes=1:N_fmodes
            M2(idx1+iter_fmodes,idx1+iter_fmodes)=-1/del_t;
            M2(idx1+iter_fmodes,idx1+iter_fmodes+2*N_fmodes)=1/del_t;
            M2(idx1+iter_fmodes+N_fmodes,idx1+iter_fmodes+N_fmodes)=-1/del_t;
            M2(idx1+iter_fmodes+N_fmodes,idx1+iter_fmodes+3*N_fmodes)=1/del_t;
        end
    else
        del_t=(ts(tt)-ts(tt-1));
        idx1=(N_t-1)*(2*N_fmodes);
        for iter_fmodes=1:N_fmodes
            M2(idx1+iter_fmodes,idx1+iter_fmodes)= 1/del_t;
            M2(idx1+iter_fmodes,idx1+iter_fmodes-(2*N_fmodes))= M2(idx1+iter_fmodes,idx1+iter_fmodes-(2*N_fmodes))-1/del_t;
            
            M2(idx1+iter_fmodes+N_fmodes,idx1+iter_fmodes+N_fmodes )=1/del_t;
            M2(idx1+iter_fmodes+N_fmodes,idx1+iter_fmodes-N_fmodes)=M2(idx1+iter_fmodes+N_fmodes,idx1+iter_fmodes-N_fmodes)-1/del_t;
            
            
        end
    end
end
 
Mbig=(M1.'*M1+1/eps.*M2.'*M2+ M3.'*M3);
Vs=Mbig\(M1.'*abs(hs).');

 
hz_slowflow=zeros(1,N_t);
slow_vars=zeros(2*N_fmodes,N_t);
for tt=1:N_t
    idx1=(tt-1)*(2*N_fmodes);
    
    hz_slowflow(tt)=mean(hs);
    slow_vars(1,tt)=mean(hs);
    for iter_fmodes=1:N_fmodes
        hz_slowflow(tt)=hz_slowflow(tt)+Vs(idx1+iter_fmodes)*abs(cos(ws(iter_fmodes)*ts(tt)))...
            + Vs(idx1+iter_fmodes+N_fmodes)*abs(sin(ws(iter_fmodes)*ts(tt)));
    slow_vars(1+iter_fmodes,tt)=Vs(idx1+iter_fmodes);        
    slow_vars(1+N_fmodes+iter_fmodes,tt)=Vs(idx1+N_fmodes+iter_fmodes);    
    end
end
 

err=sum((hz_slowflow-abs(hs)).^2);
 

end


function [hz_slowflow, slow_vars,err]= my_SF_fit_wdyn(ts, hs,slow_vars,ws,A_mat,eps,al)
N_t=ts(end);
N_fmodes=length(ws);
N_dim=(2*N_fmodes+1);
M1s=zeros(N_t,N_t,N_dim);

M3s=zeros(N_t,N_t,N_dim);


delay_dim=length(A_mat(1,:))/(N_dim);
h_vec=zeros(N_t,N_dim);




M2=spdiags([-ones(N_t,1) ones(N_t,1)],[ 0 1],N_t,N_t);
M2(end,end)=1;
M2(end,end-1)=-1;

for jj=1:N_dim
    for dd=1:N_t-delay_dim
        M3s(delay_dim+dd,dd-1+(1:delay_dim),jj)=-A_mat(jj,jj:N_dim:end);
        M3s(delay_dim+dd,dd+delay_dim,jj)=1;
    end
end


for tt=1:N_t
    if sum(ts==tt)==1
        M1s(tt,tt,1:N_dim)=1;
        h_vec(tt,:)=slow_vars(:,tt==ts);
    end
     
    
end


%hz_slowflow=zeros(1,N_t);
slow_vars=zeros(2*N_fmodes,N_t);

for jj=1:N_dim
Mbig=(M1s(:,:,jj).'*M1s(:,:,jj)+1/eps.*M2.'*M2+al.*M3s(:,:,jj).'*M3s(:,:,jj));
Vs=Mbig\(M1s(:,:,jj).'*h_vec(:,jj));

slow_vars(jj,:)=Vs;


end

tts=1:ts(end);
hz_slowflow=slow_vars(1,:);
for iter_fmodes=1:N_fmodes
hz_slowflow=hz_slowflow+slow_vars(iter_fmodes+1,:).*cos(ws(iter_fmodes).*tts)+slow_vars(N_fmodes+iter_fmodes+1,:).*sin(ws(iter_fmodes).*tts);

end
err=sum((hz_slowflow(ts)-hs).^2);
end


function [hz_slowflow, slow_vars,err]= my_SF_fit_succ(ts, hs,ws,slow_vars_IC,t0,eps)
N_t=length(ts);
N_fmodes=length(ws);

slow_vars=zeros(2*N_fmodes+1,N_t);
%slow_vars=repmat(slow_vars_IC,1,N_t);
hz_slowflow=zeros(1,N_t);
M1=zeros(1,2*N_fmodes+1);
M1(1,1)=1;
M1(1,2:N_fmodes+1)= cos(ws.*ts(1));
M1(1,N_fmodes+2:2*N_fmodes+1)= sin(ws.*ts(1));
del_t=(ts(1)-t0);

%Mbig=.*eye(2*N_fmodes+1);
Vs=1/(M1*M1.'+1/(eps*del_t^2)).*(hs(1).*M1.'+1/(eps*del_t^2).*slow_vars_IC);



hz_slowflow(1)=Vs(1);
slow_vars(1,1)=Vs(1);
for iter_fmodes=1:N_fmodes
    hz_slowflow(1)=hz_slowflow(1)+Vs(1+iter_fmodes)*cos(ws(iter_fmodes)*ts(1))...
        + Vs(1+iter_fmodes+N_fmodes)*sin(ws(iter_fmodes)*ts(1));
    slow_vars(1+iter_fmodes,1)=Vs(1+iter_fmodes);
    slow_vars(1+N_fmodes+iter_fmodes,1)=Vs(1+N_fmodes+iter_fmodes);
end
 
    
for tt=2:N_t
    M1=zeros(1,2*N_fmodes+1);
    M1(1,1)=1;
    M1(1,2:N_fmodes+1)= cos(ws.*ts(tt));
    M1(1,N_fmodes+2:2*N_fmodes+1)= sin(ws.*ts(tt));
    del_t=(ts(tt)-ts(tt-1));
  
    %Mbig=(M1*M1.'+1/(eps*del_t^2)).*eye(2*N_fmodes+1);
    Vs=1/(M1*M1.'+1/(eps*del_t^2)).*(hs(tt).*M1.'+1/(eps*del_t^2).*slow_vars(:,tt-1));

      
 
    hz_slowflow(tt)=Vs(1);
    slow_vars(1,tt)=Vs(1);
    for iter_fmodes=1:N_fmodes
        hz_slowflow(tt)=hz_slowflow(tt)+Vs(1+iter_fmodes)*cos(ws(iter_fmodes)*ts(tt))...
            + Vs(1+iter_fmodes+N_fmodes)*sin(ws(iter_fmodes)*ts(tt));
    slow_vars(1+iter_fmodes,tt)=Vs(1+iter_fmodes);        
    slow_vars(1+N_fmodes+iter_fmodes,tt)=Vs(1+N_fmodes+iter_fmodes);    
    end
 
    
end


err=sum((hz_slowflow-hs).^2);
%err2=sum((fit_wave(:,end).'-hs).^2)



 

end

function [ws_MC, err_sv]= my_Freq_MC_fit(ts,hs,err_const,N_fmodes,w0,N_MC,samplerate)

str_fmode=append('sin',num2str(N_fmodes));

 
f_low=0.5;
f_high=1.5;
if isempty(gcp('nocreate'))==1
    parpool
end
w_low_per_sampl=f_low*2*pi*samplerate;
w_high_per_sampl=f_high*2*pi*samplerate;

err_tmp=zeros(N_MC,1);
ws_tmp=zeros(N_MC+1,N_fmodes);
err_tmp(end)=err_const;
ws_tmp(end,:)=w0;


%ws=w_low_per_sampl+(w_high_per_sampl-w_low_per_sampl).*rand(N_MC,N_fmodes);
parfor iter_MC=1:N_MC
    opts=fitoptions(str_fmode);
    opts.StartPoint=zeros(1,N_fmodes*3);
    
    opts.StartPoint(2:3:end)=w_low_per_sampl+(w_high_per_sampl-w_low_per_sampl).*rand(1,N_fmodes);%[0.1 0.2 0.3 0.4];
    my_fit=fit(ts.',(hs-mean(hs)).',str_fmode,opts);
    err_tmp(iter_MC)=sum((hs-mean(hs)-my_fit(ts).').^2);
    %cost_const(jj)=err_const(jj)./err0(jj)+alpha*(2*N_fmodes)/length(hs);
    
     
     coeff2=coeffvalues(my_fit);
     ws_tmp(iter_MC,:)=sort(coeff2(2:3:end));
        
       

end

[~, idx]=min(err_tmp);
err_sv=err_tmp(idx);

ws_MC=ws_tmp(idx,:);

end

function [ws_MC, err_sv]= my_Freq_MC_fit2(ts,hs,err_const,N_fmodes,w0,N_MC,samplerate)


% From Hasselman et al. and buoy data freq sepctrum
f_low=0.02;
f_high=0.6/3;
if isempty(gcp('nocreate'))==1
    parpool
end
w_low_per_sampl=f_low*2*pi*samplerate;
w_high_per_sampl=f_high*2*pi*samplerate;

err_tmp=zeros(N_MC,1);
ws_tmp=zeros(N_MC+1,N_fmodes);
err_tmp(end)=err_const;
ws_tmp(end,:)=w0;


my_cost_fcn1=@(pars) 0;
my_cost_fcn2=@(pars) 0;
for iter_fmodes=1:N_fmodes
    idx1=3*(iter_fmodes-1);
my_cost_fcn1=@(pars) my_cost_fcn1(pars(1:idx1))+pars(idx1+1).*cos(pars(idx1+3).*ts)+pars(idx1+2).*sin(pars(idx1+3).*ts);

my_cost_fcn2=@(pars) my_cost_fcn2(pars(1:idx1))-pars(idx1+1)*pars(idx1+3).*sin(pars(idx1+3).*ts)+pars(idx1+2)*pars(idx1+3).*cos(pars(idx1+3).*ts);

end
my_cost_fcn=@(pars) [my_cost_fcn1(pars)-hs my_cost_fcn2(pars)];
opts = optimoptions('lsqnonlin','MaxFunctionEvaluations',10^5,'FunctionTolerance',10^-6,'MaxIterations',10^4,'Display','Off');%


%ws=w_low_per_sampl+(w_high_per_sampl-w_low_per_sampl).*rand(N_MC,N_fmodes);
parfor iter_MC=1:N_MC
    %     opts=fitoptions(str_fmode);
    %     opts.StartPoint=zeros(1,N_fmodes*3);
    %
    ws_start=w_low_per_sampl+(w_high_per_sampl-w_low_per_sampl).*rand(1,N_fmodes);%[0.1 0.2 0.3 0.4];
    %     my_fit=fit(ts.',(hs-mean(hs)).',str_fmode,opts);
    pars0=zeros(3*N_fmodes,1);
    pars0(3:3:end)=ws_start;
    [coeffs, err_tmp(iter_MC)]=lsqnonlin(my_cost_fcn,pars0,[],[],opts)
    %err_tmp(iter_MC)=sum((hs-mean(hs)-my_fit(ts).').^2);
    %cost_const(jj)=err_const(jj)./err0(jj)+alpha*(2*N_fmodes)/length(hs);
    
    
    ws_tmp(iter_MC,:)=sort(coeffs(3:3:end));
    
    
    
end

[~, idx]=min(err_tmp);
err_sv=err_tmp(idx);

ws_MC=ws_tmp(idx,:);

end


function hz_slowflow_full=SF2fulltime(t_full,ts,slow_vars,ws,N_fmodes)
      hz_slowflow_full=zeros(1,length(t_full));
      for tt=1:length(t_full)
          
          %    hz_slowflow_tmp(tt)=sum(tmp_coeffs(:,1));%
          hz_slowflow_full(tt)=interp1(ts,slow_vars(1,:),t_full(tt),'linear','extrap');
          for iter_fmodes=1:N_fmodes
              hz_slowflow_full(tt)=hz_slowflow_full(tt)+interp1(ts,slow_vars(1+iter_fmodes,:),t_full(tt),'linear','extrap')*cos(ws(iter_fmodes)*t_full(tt))...
                  +  interp1(ts,slow_vars(1+N_fmodes+iter_fmodes,:),t_full(tt),'linear','extrap')*sin(ws(iter_fmodes)*t_full(tt));
          end
      end 
end

function hz_slowflow_full=SF2fulltime2(t_full,ts,slow_vars,ws,N_fmodes)
      hz_slowflow_full=zeros(1,length(t_full));
      for tt=1:length(t_full)
          
          %    hz_slowflow_tmp(tt)=sum(tmp_coeffs(:,1));%
          hz_slowflow_full(tt)=interp1(ts,slow_vars(1,:),t_full(tt),'linear','extrap');
          for iter_fmodes=1:N_fmodes
              hz_slowflow_full(tt)=hz_slowflow_full(tt)+interp1(ts,slow_vars(1+iter_fmodes,:),t_full(tt),'linear','extrap')*abs(cos(ws(iter_fmodes)*t_full(tt)))...
                  +  interp1(ts,slow_vars(1+N_fmodes+iter_fmodes,:),t_full(tt),'linear','extrap')*abs(sin(ws(iter_fmodes)*t_full(tt)));
          end
      end 
end

function [h_cmp, t_cmp]=my_wave_height_compare(true_wave_height,true_wave_time,wave_height,wave_time)

h_cmp=zeros(1,length(true_wave_height));
t_cmp=zeros(1,length(true_wave_height));

for iter_waves=1:length(true_wave_time)
    
    [~,idx]=min(vecnorm([true_wave_height(iter_waves);true_wave_time(iter_waves)]-[wave_height;wave_time]));
    
    h_cmp(iter_waves)=wave_height(idx);
    t_cmp(iter_waves)=wave_time(idx);
 
end
    

 

end



%   