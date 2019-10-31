%% 30CST
% Written by Lindsey Tulipani - lindsey.tulipani@uvm.edu
% October 31, 2019
% This code computes accelerometer based metrics (including number of stands 
% performed) using thigh and chest sensors to analyze the 30 sec repeated 
% sit-to-stand functional test. 
% *Note: This code does not take into consideration false starts/false stands. 

% Load example data
load('time_s.mat')
load('time_cs.mat')
load('th_acc_stat.mat')
load('th_acc_cs.mat')
load('ch_acc_stat.mat')
load('ch_acc_cs.mat')

%Rotate sensor frame to thigh frame
accel_stat = mean([th_acc_stat(:,1),th_acc_stat(:,2),th_acc_stat(:,3)]',2);
th_acc_stat_n = accel_stat/vecnorm(accel_stat);
z = [0 0 1]'; 
[th_R] = getrot(th_acc_stat_n,z,'dcm');

thigh_accel = th_R*[th_acc_cs(:,1),th_acc_cs(:,2),th_acc_cs(:,3)]';
thigh_accel = thigh_accel';

%Rotate sensor frame to chest frame
g = mean([ch_acc_stat(:,1),ch_acc_stat(:,2),ch_acc_stat(:,3)]',2);
g=g';
z = g ./ norm(g);
x = cross([0,0,1],z)./norm(cross([0,0,1],z));
y = cross(z,x);
Cr = [x; y; z];
chest_accel = Cr*[ch_acc_cs(:,1),ch_acc_cs(:,2),ch_acc_cs(:,3)]';
chest_accel = [chest_accel]';

%Compute frequency characteristics of chest 
Fs = 1/mean(diff(time_cs));
[p,F] = pwelch(chest_accel(:,2),rectwin(round(5*Fs)),[],4096,Fs);
freq = [F,p];

[freq_pks,locs_f_pks] = findpeaks(p);
freqs = F(locs_f_pks(F(locs_f_pks)>=0.4));
[~,imax] = max(freq_pks(F(locs_f_pks)>=0.4));
ch_f=freqs(imax);

% Applylow-pass butterworth filter with a cutoff Hz 
n = 3; % filter order 
Wp = ch_f/(Fs/2); % defining cutoff frequency in rad/samp
[b,a] = butter(n,Wp,'low'); % 
ch_yaccel = filtfilt(b,a,chest_accel(:,2)); 
ch_zaccel = filtfilt(b,a,chest_accel(:,3)); 

%Compute chest angle
ch_theta = rad2deg(atan2(ch_yaccel,ch_zaccel));

%Compute frequency characteristics of thigh
[p_th,F_th] = pwelch(thigh_accel(:,3),rectwin(round(3*Fs)),[],4096,Fs);
[tfreq_pks,tlocs_f_pks] = findpeaks(p_th);
th_f = F_th(tlocs_f_pks(2));

% Applylow-pass butterworth filter with a cutoff Hz
n = 3; % filter order 
Wp = th_f/(Fs/2); 
[b,a] = butter(n,Wp,'low'); % 
th_zaccel = filtfilt(b,a,thigh_accel(:,3)); 

%Location of stands
[pks_stand,locs_stand] = findpeaks(th_zaccel);

%Number of repetitions performed
num_reps = length(pks_stand);

%Location of sits
locs_sit=zeros(1,num_reps-1);
for k=1:(num_reps-1)
    [pks,locs] = findpeaks(-th_zaccel);
        if locs(1) < locs_stand(1)
            locs_sit(k) = locs(k+1);
            pks_sit(k) = pks(k+1);
        else locs_sit(k) = locs(k);
            pks_sit(k) = pks(k);
        end
end

%Compute sit to stand transition times
sit_stand_time = zeros(1,num_reps);

for k=1:num_reps
if k==1 
    sit_stand_time(k) = time_cs(locs_stand(k))-time_cs(k);
end
if k == num_reps
    sit_stand_time(k) = time_cs(locs_stand(k))-time_cs(locs_sit(k-1));
end

if 1 < k && k < num_reps
    sit_stand_time(k) = time_cs(locs_stand(k))-time_cs(locs_sit(k-1));
end
end

%Calc avg/sd/cv si-st time
avg_sit_stand_t = mean(sit_stand_time);
std_sit_stand_t = std(sit_stand_time);
cov_sit_stand_t = std_sit_stand_t / avg_sit_stand_t;

%Compute stand to sit transition times
stand_sit_time = zeros(1,num_reps-1);
for k=1:(num_reps-1)
    stand_sit_time(k) = time_cs(locs_sit(k))-time_cs(locs_stand(k));
end
%Calc avg/sd/cv of st-si times
avg_stand_sit_t = mean(stand_sit_time);
std_stand_sit_t = std(stand_sit_time);
cov_stand_sit_t = std_stand_sit_t / avg_stand_sit_t;

%Compute Trunk Flex during sit to stand

trflex_si_st = zeros(1,num_reps);
for k=1:num_reps
    if k==1 
    trflex_si_st(k) = min(ch_theta(1:locs_stand(1)));
    end
    if 1 < k && k <= num_reps
    trflex_si_st(k) = min(ch_theta(locs_sit(k-1):locs_stand(k)));
    end
end
%Calc avg/sd/cv of si-st trunk flex
avg_trflex_si_st = mean(trflex_si_st);
std_trflex_si_st = std(trflex_si_st);
cov_trflex_si_st = std_trflex_si_st / avg_trflex_si_st;

%Band pass filter accel data

% Apply high-pass butterworth filter thigh
n = 3; % filter order 
Wp = 5/(Fs/2); % defining cutoff frequency of 5 Hz in rad/samp
[b,a] = butter(n,Wp,'high'); 
th_xaccelHP = filtfilt(b,a,thigh_accel(:,1)); 
th_yaccelHP = filtfilt(b,a,thigh_accel(:,2)); 
th_zaccelHP = filtfilt(b,a,thigh_accel(:,3)); 

n = 3; % filter order 
Wp = 30/(Fs/2); % defining cutoff frequency of 30 Hz in rad/samp
[b,a] = butter(n,Wp,'low'); 
th_xaccelHP = filtfilt(b,a,th_xaccelHP); 
th_yaccelHP = filtfilt(b,a,th_yaccelHP); 
th_zaccelHP = filtfilt(b,a,th_zaccelHP);

%Find inflection points to identify ranges of interest
diff_az = diff(th_zaccel); %find deriv of accel_z 
[pks_tstand,locs_tstand] = findpeaks(diff_az); %peak represents inflec point sit to stand

%find inflection points from stand to sit
[pks_tsit,locs_tsit] = findpeaks(-diff_az); %peak represents inflec point st-si

%Find peak thigh accel (magnitude) of the si-st transition

pk_si_st_vacc = zeros(1,num_reps);
pk_si_st_APacc = zeros(1,num_reps);
pk_si_st_MLacc = zeros(1,num_reps);

abs_th_xaccelHP = abs(th_xaccelHP); %absolute value of HP filtered thigh signal
abs_th_yaccelHP = abs(th_yaccelHP); 
abs_th_zaccelHP = abs(th_zaccelHP); 

for k=1:num_reps
 if k==1
    pk_si_st_vacc(k) = max(abs_th_zaccelHP(1:locs_tstand(k),1));
    pk_si_st_APacc(k) = max(abs_th_yaccelHP(1:locs_tstand(k),1));
    pk_si_st_MLacc(k) = max(abs_th_xaccelHP(1:locs_tstand(k),1));
 end
if 1 < k && k <= (num_reps)
    pk_si_st_vacc(k) = max(abs_th_zaccelHP(locs_sit(k-1):locs_tstand(k),1));
    pk_si_st_APacc(k) = max(abs_th_yaccelHP(locs_sit(k-1):locs_tstand(k),1));
    pk_si_st_MLacc(k) = max(abs_th_xaccelHP(locs_sit(k-1):locs_tstand(k),1));
end
end

pk_st_si_vacc = zeros(1,num_reps-1);
pk_st_si_APacc = zeros(1,num_reps-1);
pk_st_si_MLacc = zeros(1,num_reps-1);
trflex_st_si = zeros(1,num_reps-1);

%Find peak thigh accel (magnitude) of the st-si transition
for k=1:(num_reps-1)
    pk_st_si_vacc(k) = max(abs_th_zaccelHP(locs_tsit(k):locs_sit(k),1));
    pk_st_si_APacc(k) = max(abs_th_yaccelHP(locs_tsit(k):locs_sit(k),1));
    pk_st_si_MLacc(k) = max(abs_th_xaccelHP(locs_tsit(k):locs_sit(k),1));
    trflex_st_si(k) = max(ch_theta(locs_tsit(k):locs_sit(k),1));
end

%Calc avg/sd/cv of st-si trunk flex
avg_trflex_st_si = mean(trflex_st_si);
std_trflex_st_si = std(trflex_st_si);
cov_trflex_st_si = std_trflex_st_si / avg_trflex_st_si;

%Avg, Std and CV of metrics

%Si-st
avg_peak_vacc_si_st = mean(pk_si_st_vacc);
std_peak_vacc_si_st = std(pk_si_st_vacc);
cov_peak_vacc_si_st =  std_peak_vacc_si_st / avg_peak_vacc_si_st;

avg_peak_MLacc_si_st = mean(pk_si_st_MLacc);
std_peak_MLacc_si_st = std(pk_si_st_MLacc);
cov_peak_MLacc_si_st =  std_peak_MLacc_si_st / avg_peak_MLacc_si_st;

avg_peak_APacc_si_st = mean(pk_si_st_APacc);
std_peak_APacc_si_st = std(pk_si_st_APacc);
cov_peak_APacc_si_st =  std_peak_APacc_si_st / avg_peak_APacc_si_st;

%St-si
avg_peak_vacc_st_si = mean(pk_st_si_vacc);
std_peak_vacc_st_si = std(pk_st_si_vacc);
cov_peak_vacc_st_si =  std_peak_vacc_st_si/ avg_peak_vacc_st_si;

avg_peak_MLacc_st_si = mean(pk_st_si_MLacc);
std_peak_MLacc_st_si = std(pk_st_si_MLacc);
cov_peak_MLacc_st_si =  std_peak_MLacc_st_si/ avg_peak_MLacc_st_si;

avg_peak_APacc_st_si = mean(pk_st_si_APacc);
std_peak_APacc_st_si = std(pk_st_si_APacc);
cov_peak_APacc_st_si =  std_peak_APacc_st_si/ avg_peak_APacc_st_si;

% Define table to report results
h1 = {'num_reps','avg_trflex_si_st','cov_trflex_si_st','avg_trflex_st_si','cov_trflex_st_si','avg_sit_stand_t','cov_sit_stand_t',...
    'avg_stand_sit_t','cov_stand_sit_t','avg_peak_vacc_si_st','cov_peak_vacc_si_st','avg_peak_vacc_st_si','cov_peak_vacc_st_si',...
    'avg_peak_MLacc_si_st', 'cov_peak_MLacc_si_st', 'avg_peak_MLacc_st_si','cov_peak_MLacc_st_si', 'avg_peak_APacc_si_st', 'cov_peak_APacc_si_st'};

mets = {num_reps,avg_trflex_si_st,cov_trflex_si_st,avg_trflex_st_si,cov_trflex_st_si,avg_sit_stand_t,cov_sit_stand_t,...
    avg_stand_sit_t,cov_stand_sit_t,avg_peak_vacc_si_st,cov_peak_vacc_si_st,avg_peak_vacc_st_si,cov_peak_vacc_st_si,...
    avg_peak_MLacc_si_st, cov_peak_MLacc_si_st, avg_peak_MLacc_st_si,cov_peak_MLacc_st_si, avg_peak_APacc_si_st, cov_peak_APacc_si_st};

T1 = cell2table(mets,'VariableNames',h1);

% Write table to command window
disp(T1);


%% Supporting Functions
function [ r ] = getrot( v1, v2, type)
%Reed Gurchiek, 2017
%   getrot finds the rotation operator of type 'type' which takes v 
%   measured in frame 1 (v1) and expresses it in frame 2 (v2) if v2 is a 
%   3-dimensional vector.  Otherwise, it constructs an angle-axis
%   rotator where v1 is the axis and v2 is the angle if v2 is
%   1-dimensional.  In this case, consider the axis (in frame 1) and angle
%   that one would use to rotate frame 1 to align with frame 2
%
%-----------------------------INPUTS---------------------------------------
%
%   v1, v2:
%       vectors 1 and 2. 3xn matrix of column vectors. v1 is v measured in
%       frame 1 and v2 is v measured in frame 2.
%       OR
%       rotation axis (v1: 3xn matrix of column vectors with unit norm) and
%       rotation angle (v2: 1xn array of rotation angles)
%
%   type:
%       string specifying type of rotation operator.  Either 'dcm' for
%       direction cosine matrix or 'q' for quaternion.
%
%----------------------------OUTPUTS---------------------------------------
%
%   r:
%       rotation operator which takes v1 to v2 of type 'type' or described
%       by the axis-angle combo v1 & v2.
%
%--------------------------------------------------------------------------

%% getrot

%verify proper inputs
[v1r,v1c] = size(v1);
[v2r,v2c] = size(v2);
if v1r ~= 3
    error('v1 must be 3xn matrix of column vectors')
elseif v1c ~= v2c
    error('v1 and v2 must have same number of columns')
elseif v2r ~= 1 && v2r ~= 3
    error('v2 must either be 3xn (if a vector) or 1xn (if an angle)')
end

%if v2 is a vector
if v2r == 3
    
    %get axis of rotation
    axis = cross(v2,v1)./vecnorm(cross(v2,v1));
    
    %get angle
    angle = acos(dot(v2,v1)./(vecnorm(v2).*vecnorm(v1)));

%if v2 is 1D array of angles
elseif v2r == 1
    
    %axis and angle given
    axis = v1;
    angle = v2;
    
end

%if quaternion
if strcmpi(type,'q')
    
    %construct quaternion
    r = [repmat(sin(angle/2), [3 1]).*axis; cos(angle/2)];
    
%if dcm
elseif strcmpi(type,'dcm')
    
    %construct dcm (euler formula: R(n,a) = I - s(a)*[nx] + (1-c(a))*[nx]^2)
    r = zeros(3,3,v1c);
    I3 = eye(3);
    for k = 1:v1c
        r(:,:,k) = I3 - sin(angle(k))*skew(axis(:,k)) + (1-cos(angle(k)))*skew(axis(:,k))*skew(axis(:,k));
    end
end

end

function [ Vx ] = skew( V )
%Reed Gurchiek, 2017
%   skew takes a 3xn matrix of column vectors and returns a 3x3xn skew
%   symmetric matrix for each column vector in V such that Vx(3,3,i)*p =
%   cross(V(:,i),p).
%
%---------------------------------INPUTS-----------------------------------
%
%   V:
%       3xn matrix of column vectors.
%
%--------------------------------OUTPUTS-----------------------------------
%
%   Vx:
%       3x3xn skew symmetric matrices.
%
%--------------------------------------------------------------------------

%% skew

%verify proper inputs
[Vr,Vc] = size(V);
if Vc == 3 && Vr ~= 3
    V = V';
elseif Vr ~= 3 && Vc ~= 3
    error('V must have 3 rows or 3 columns')
end

%for each vector
[~,n] = size(V);
Vx = zeros(3,3,n);
for k = 1:n
    
    %get skew
    Vx(:,:,k) = [   0    -V(3,k)  V(2,k);...
                  V(3,k)    0    -V(1,k);...
                 -V(2,k) V(1,k)     0   ];
             
end



end

