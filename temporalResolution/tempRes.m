%% create spatial map
% circular object
disp 'Generating phantom'
dim = 50; 
fov = 20;       % cm
rdisc1 = 0.7*fov/2;    % disc1 radius is % of FOV/2
rdisc2 = 0.9*fov/2;    % disc2 radius is % of FOV/2

x = -fov/2:fov/dim:fov/2-fov/dim;
y = -fov/2:fov/dim:fov/2-fov/dim;

[X,Y] = meshgrid(x,y);

mag = zeros(dim);
mag(X.^2 + Y.^2 <= rdisc1^2) = 0.5;
mag(X.^2 + Y.^2 <= rdisc2^2) = mag(X.^2 + Y.^2 <= rdisc2^2) + 0.5;

    figure; 
    % plot map
    subplot(2,1,1); 
    imagesc(x,y,mag); axis image; colorbar; %colormap(gray);
    title('Spatial map','FontSize',14); 
    xlabel('x [cm]','FontSize',14); ylabel('y [cm]','FontSize',14); 
    % plot profile
    subplot(2,1,2); 
    plot(x,mag(:,length(mag)/2),'Linewidth',2); 
    title('Profile along x','FontSize',14); 
    xlabel('x [cm]','FontSize',14); ylabel('Intensity','FontSize',14);

%% create a pulse
s = 0.3; % width
Z = exp( -(X.^2+Y.^2)/(2*s^2) ); 
Z0 = Z/max(Z(:)); % normalized pressure profile


    figure;
    imagesc(x,y,Z0); axis image; colorbar;
    title('Pressure pulse','FontSize',14); 
    xlabel('x [cm]','FontSize',14); ylabel('y [cm]','FontSize',14);

%% temperature change simulation

% specify vector of timepoints
tmax = 90;              % duration of temperature simulation [s] 
tstep = 20e-3;          % step size [s]
Ntp = tmax/tstep;       % number of timepoints
tpulse = 25;            % duration of pulse [s]

t = 0:tmax/(Ntp-1):tmax;
D = 1;                  % diffusivity           %%%%% UPDATE THIS
alpha = 6;           % chemical shift coefficient [ppm/degreeC] 
                        % (see Todd et al. 2009)
gamma = 42.57e6*2*pi;   % gyromagnetic ratio [rad Hz/T]
B0 = 3;                 % field strenght [T]
TE = 12e-3;             % echo time [s]

Tmap = 37 + zeros(size(Z0,1),size(Z0,2),length(t));  % body initial temperature

% Model predictive filtering (Todd, Payne, Parker 2010)
% Tn+1 = Tn + (k/rho/C*del2Tn - W/rhoTn + un/rho/C*Q)*dt
% k: thermal conductivity (W/m/deg C)
% rho: tissue density (kg/m3)
% C: specific heat (J/kg/deg C)
% W = perfusion; ignore
% un: input power (watts)
rho = 1000; % kg/m3, tissue density
C = 4186; % J/kg/deg C, specific heat
u = 21500000; % Watts - why so huge??
k = 0.4*100; % Watts/m/deg C - why so huge?? 0.4 is from Todd

% hypothetical temperature distribution
% estimate change over space and time
disp 'Simulating bioheat equation'
for idx = 2:length(t)
    if t(idx) <= tpulse
        deltatemp = (k*del2(Tmap(:,:,idx-1),fov/100/dim,fov/100/dim)+u*Z0)*tstep/(rho*C);
        Tmap(:,:,idx) = Tmap(:,:,idx-1) + deltatemp;
    else
        deltatemp = k*del2(Tmap(:,:,idx-1),fov/100/dim,fov/100/dim)*tstep/(rho*C);
        Tmap(:,:,idx) = Tmap(:,:,idx-1) + deltatemp;
    end
end


figure; h = plot(t,squeeze(Tmap(dim/2+1,dim/2+1,:)),'linewidth',2); 
title('Peak Temperature vs Time','FontSize',14); 
xlabel('Time [s]','FontSize',14); ylabel('Temperature [deg C]','FontSize',14);

return

% truncate it to +/- 10 seconds around the peak
[~,maxi] = max(sqz(Tmap(dim/2+1,dim/2+1,:)));
peakInds = (t > (t(maxi)-15)) & (t < (t(maxi)+15));
Tmap = Tmap(:,:,peakInds);
nrow = 4;

    
% phase map
phi = (Tmap-37)*0.01*10^-6*gamma*B0*TE;
 
cmax = max(max(phi(:)));
clim = [0 cmax];

% create complex images
cplxdata = repmat(mag,[1 1 size(Tmap,3)]).*exp(1i*phi);

%% direct temperature sampling
disp 'Direct Temp Sampling'

% sample the temp curve with different temporal resolutions; 
% record highest error among grid shifts
tempCurve = squeeze(Tmap(dim/2+1,dim/2+1,:));
res = 10:10:500; % integer steps
worstErrSingle = zeros(size(res));
worstErrAve = zeros(size(res));
bestErrSingle = Inf + zeros(size(res)); % trivial; should be zero
bestErrAve = Inf + zeros(size(res)); % synchronized error
%[~,iMax] = max(tempCurve(:));
for ii = 1:length(res)
    % loop over offsets, record highest error for both window-averaged and
    % single-point
    for jj = 1:res(ii) 
        sampledSingle = tempCurve(jj:res(ii):end);
        maxErrSingle = max(tempCurve)-max(sampledSingle);
        if maxErrSingle > worstErrSingle(ii)
            worstErrSingle(ii) = maxErrSingle;
        end
        if maxErrSingle < bestErrSingle(ii)
            bestErrSingle(ii) = maxErrSingle;
        end
        % average every consecutive res(ii) time points to get sampledAve
        sampledAve = tempCurve(jj:end);
        sampledAve = sampledAve(1:floor(length(sampledAve)/res(ii))*res(ii));
        sampledAve = mean(reshape(sampledAve,[res(ii) length(sampledAve)/res(ii)]),1);
        maxErrAve = max(tempCurve) - max(sampledAve);
        if maxErrAve > worstErrAve(ii)
            worstErrAve(ii) = maxErrAve;
        end
        if maxErrAve < bestErrAve(ii)
            bestErrAve(ii) = maxErrAve;
        end
    end
    
end

figure;hold on
plot(res*tstep,worstErrSingle)
plot(res*tstep,worstErrAve)
plot(res*tstep,bestErrAve)
ylabel 'Degrees Celsius'
xlabel 'Seconds'
legend('Single-point error, worst case','Averaged error, worst case','Averaged error, synchronized');

return

%% 2DFT simulation
disp '2DFT Temp Sampling'

% loop over TRs/resolutions and offsets, calculate worst-case and
% synchronized error
res2DFT = 1:10; % integer time steps between each k-space line (corresponds to 
% frame rates of 1 - 10 seconds)
bestErr2DFT = Inf + zeros(size(res2DFT));
worstErr2DFT = zeros(size(res2DFT));
tempCurve = squeeze(Tmap(dim/2+1,dim/2+1,:)); % max temp curve
% ft all the images
for ii = 1:size(cplxdata,3);
    cplxkData(:,:,ii) = ft2(cplxdata(:,:,ii));
end
for ii = 1:length(res2DFT) % loop over TRs/frame rates
    ii
    for jj = 1:dim*res2DFT(ii) % loop over temporal offsets for each frame rate
        dataAll = cplxkData(:,:,jj:res2DFT(ii):end); % images at time points we will extract k-data from
        dataAll = dataAll(:,:,1:floor(size(dataAll,3)/dim)*dim); % trim to size
        % grab one k-space line from each time point
        dataSampled = zeros([dim dim size(dataAll,3)/dim]);
        for kk = 1:dim
            dataSampled(:,kk,:) = dataAll(:,kk,kk:dim:end);
        end
        % ift dataSampled
        for kk = 1:size(dataSampled,3)
            dataSampled(:,:,kk) = ift2(dataSampled(:,:,kk));
        end
        % get the temp in the middle
        sampledTemp = angle(squeeze(dataSampled(dim/2+1,dim/2+1,:)));
        sampledTemp = sampledTemp./(0.01*10^-6*gamma*B0*TE) + 37;
        % get difference between true and measured max temps
        maxErr2DFT = max(tempCurve) - max(sampledTemp);
        %[ii jj]
        %figure(100);plot(sampledTemp);pause
        if maxErr2DFT > worstErr2DFT(ii)
            worstErr2DFT(ii) = maxErr2DFT;
        end
        if maxErr2DFT < bestErr2DFT(ii)
            bestErr2DFT(ii) = maxErr2DFT;
        end
    end
end
figure;hold on
plot(res2DFT*dim*tstep,worstErr2DFT)
plot(res2DFT*dim*tstep,bestErr2DFT)
ylabel 'Degrees Celsius'
xlabel 'Seconds'
legend('2DFT error, worst case','2DFT error, synchronized');

return
%% spiral simulation
if ~exist('Gmri') || ~exist('qpwls_pcg')
    error 'Need Jeff Fesslers Image Recon Toolbox installed. Download at http://web.eecs.umich.edu/~fessler/code/index.html'
end
if ~exist('vds')
    error 'Need variable density spiral function from Brian Hargreaves. Download at http://mrsrl.stanford.edu/~brian/vdspiral/'
end
disp 'Spiral Temp Sampling'

% get a dim-shot spiral trajectory
addpath vdspiral/
smax = 15000;	 % 150 T/m/s
gmax = 4;	 % G/cm
dt = 4e-6; % dwell time of spiral in seconds
N = dim; % Interleaves
Fcoeff = [fov 0]; 	% FOV starts at 20 cm; doesn't decrease.
res = fov/dim/sqrt(2); % cm, resolution
rmax = 1/2/res;	% cm^(-1)
[k,g,s,time,r,theta] = vds(smax,gmax,dt,N,Fcoeff,rmax);
kAll = k(:)*exp(1i*2*pi/dim*(0:dim-1)); % duplicate to all shots

% build a system matrix
mask = true(dim);
G = Gmri([real(kAll(:)) imag(kAll(:))],mask,'fov',fov);

% loop over TRs/resolutions and offsets, calculate worst-case and
% synchronized error
resSpiral = 1:10; % integer time steps between each k-space line (corresponds to 
% 1 - 10 seconds)
bestErrSpiral = Inf + zeros(size(resSpiral));
worstErrSpiral = zeros(size(resSpiral));
tempCurve = squeeze(Tmap(dim/2+1,dim/2+1,:)); % max temp curve
% ft all the images
spiralData = zeros(length(k),dim,size(cplxdata,3));
parfor ii = 1:size(cplxdata,3);
    spiralData(:,:,ii) = reshape(G*col(cplxdata(:,:,ii)),[length(k) dim]);
end
nSpiralIters = 15; % spiral CG image recon iterations
for ii = 1:length(resSpiral) % loop over TRs/frame rates
    ii
    for jj = 1:dim*resSpiral(ii) % loop over temporal offsets for each frame rate
        dataAll = spiralData(:,:,jj:resSpiral(ii):end); % images at time points we will extract k-data from
        dataAll = dataAll(:,:,1:floor(size(dataAll,3)/dim)*dim); % trim to size
        % grab one spiral shot from each time point
        dataSampled = zeros([length(k) dim size(dataAll,3)/dim]);
        for kk = 1:dim
            dataSampled(:,kk,:) = dataAll(:,kk,kk:dim:end);
        end
        % recon each dynamic 
        dataRecond = zeros([dim dim size(dataSampled,3)]);
        parfor kk = 1:size(dataSampled,3)
            [xS,info] = qpwls_pcg(zeros(dim*dim,1),G,1,col(dataSampled(:,:,kk)),0,0,1,nSpiralIters,true(dim));
            dataRecond(:,:,kk) = reshape(xS(:,end),[dim dim]);
        end
%         dataRecond = zeros([dim dim size(dataAll,3)/dim]);
%         for kk = 1:size(dataSampled,3)
%             [xS,info] = qpwls_pcg(zeros(dim*dim,1),G,1,col(dataAll(:,:,1+(kk-1)*dim)),0,0,1,nSpiralIters,true(dim));
%             dataRecond(:,:,kk) = reshape(xS(:,end),[dim dim]);
%         end
        % get the temp in the middle
        sampledTemp = angle(squeeze(dataRecond(dim/2+1,dim/2+1,:)));
        sampledTemp = sampledTemp./(0.01*10^-6*gamma*B0*TE) + 37;
        % get difference between true and measured max temps
        maxErrSpiral = max(tempCurve) - max(sampledTemp);
        %[ii jj]
        %figure(100);plot(sampledTemp);pause
        if maxErrSpiral > worstErrSpiral(ii)
            worstErrSpiral(ii) = maxErrSpiral;
        end
        if maxErrSpiral < bestErrSpiral(ii)
            bestErrSpiral(ii) = maxErrSpiral;
        end
    end
end
figure;hold on
plot(resSpiral*dim*tstep,worstErrSpiral)
plot(resSpiral*dim*tstep,bestErrSpiral)
ylabel 'Degrees Celsius'
xlabel 'Seconds'
legend('Spiral error, worst case','Spiral error, synchronized');


%% plot 2DFT and spiral together
figure;hold on
plot(res2DFT*dim*tstep,worstErr2DFT)
plot(res2DFT*dim*tstep,bestErr2DFT)
plot(resSpiral*dim*tstep,worstErrSpiral)
plot(resSpiral*dim*tstep,bestErrSpiral)
plot(1:10,ones(10,1),'k--');
ylabel 'Peak Error (Degrees Celsius)'
xlabel 'Frame Rate (Seconds)'
legend('2DFT error, worst case','2DFT error, synchronized',...
    'Spiral error, worst case','Spiral error, synchronized');
axis([1 10 0 5]);
axis square


