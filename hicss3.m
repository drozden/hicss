% OUTPUTS = HICSS3(NUMTRIALS,NUMMOTES,FILENAME,FLAG_PLOT,FLAG_DISP)
function [outputs] = hicss3(varargin)
if nargin < 1
    numtrials = 10000;
else
    numtrials = varargin{1};
end
if nargin < 2
    minnummotes = 20;
    maxnummotes = 50;
else
    maxnummotes = varargin{2};
    if length(maxnummotes) > 1
        minnummotes = maxnummotes(1);
        maxnummotes = maxnummotes(2);
    end 
end
if nargin < 3
    filename = 'area3.JPG';
else
    filename = varargin{3};
end
if nargin < 4
    flag_plot = false;
else
    flag_plot = varargin{4};
end
if nargin < 5
    flag_disp = false;
else
    flag_disp = varargin{5};
end
c = physconst('lightspeed');
freq = 2.4e9;
power_tx = 8; %dBm
sensitivity_rx = -103; %dBm
loss_wall = 5; %dB
% 200 ft ~= 114 pixels
pix2ft = 200/114; ft2m = 0.3048;
sensor_range = 100/pix2ft; % 100 ft
% Fog Server Locations
% TechX Lab ~= (1031, 236)
% CoB Office ~= (1324, 330)
fogsrv_col = [1031, 1324];
fogsrv_row = [ 236,  330];
numfogsrv = length(fogsrv_col);

[I2,nr,nc] = fix_gmap(filename,fogsrv_col,fogsrv_row,flag_plot);
if flag_plot, figure(1); end
outputs(numtrials) = struct('mote_row',0,'mote_col',0,'made_connection',...
    false,'num_hops',0,'next_hop',0);
fprintf('Starting %d trials for %2d to %2d motes.\n', ...
    numtrials,minnummotes,maxnummotes);
for o = minnummotes:maxnummotes
nummotes = o;
for n = 1:numtrials
mote_row = zeros(nummotes,1); mote_col = mote_row;
dist_m = zeros([nummotes,numfogsrv]); apathloss = dist_m; 
power_rx = zeros([nummotes,numfogsrv]); 
diff_col = zeros(numfogsrv,1); diff_row = diff_col; npathpts = diff_col;
made_connection = false(nummotes,1); make_connection = false(nummotes,1);
num_hops = zeros(nummotes,1); next_hop = zeros(nummotes,1);
path_x = cell(numfogsrv,1); path_y = path_x;
circle_col = zeros(360,1); circle_row = circle_col; 
I3 = false(nr,nc);
for k = 1:nummotes 
    % For each trial, generate the random locations of each sensor mote 
    mote_row(k) = randi(nr);
    mote_col(k) = randi(nc);
    while (I2(mote_row(k),mote_col(k))==0) % Ensure each mote is outdoors
        mote_row(k) = randi(nr);
        mote_col(k) = randi(nc);
    end
    % Find (& Show) range of accurate sensor measurement
    a = linspace(0,2*pi,360);
    for p = 1:length(a)
        circle_col(p) = min(nc,max(1,round(mote_col(k)+sensor_range*cos(a(p)))));
        circle_row(p) = min(nr,max(1,round(mote_row(k) + sensor_range * sin(a(p)))));
    end
    circle_pts = unique([circle_col,circle_row],'rows');
    if flag_plot, plot(circle_pts(:,1),circle_pts(:,2),'k.'); end
    for q = min(circle_col):max(circle_col)
       r = circle_row(circle_col == q);
       I3(min(r):max(r),q) = true;
    end
    for m = 1:numfogsrv
        % Estimate distance between random mote location & each fog server
        diff_col(m) = abs(mote_col(k)-fogsrv_col(m));
        diff_row(m) = abs(mote_row(k)-fogsrv_row(m));
        dist_m(k,m) = pix2ft*ft2m*sqrt(diff_col(m)^2+diff_row(m)^2);
        % Estimate the propagation loss at this distance
        apathloss(k,m) = fspl(dist_m(k,m),c./freq);
        % Estimate the max possible received power at the fog server
        power_rx(k,m) = power_tx - apathloss(k,m);
        % Find straight path b/w mote & fog unit, estimate #obstructions
        npathpts(m) = max(diff_col(m),diff_row(m))+1;
        fit_x = [mote_col(k),fogsrv_col(m)]; 
        fit_y = [mote_row(k),fogsrv_row(m)];
        pc = polyfit(fit_x,fit_y,1);
        if (diff_col(m) >= diff_row(m))
            path_x{m} = min(fit_x):max(fit_x);
            path_y{m} = round(polyval(pc,path_x{m}));
        else
            if (min(fit_x) == max(fit_x))
                path_x{m} = fit_x(1)+zeros(1,npathpts(m));
                path_y{m} = min(fit_y):max(fit_y);
            else
            path_x{m} = linspace(min(fit_x),max(fit_x),npathpts(m));
            path_y{m} = round(polyval(pc,path_x{m}));
            assert(min(path_y{m})==min(fit_y),'1st y-pt in path is not min');
            assert(max(path_y{m})==max(fit_y),'Last y-pt in path is not max');
            path_x{m} = round(path_x{m});
            end
        end
        for p = 1:npathpts(m)
            % Decrease Rx'd power for each impediment in path to fog server
            if (I2(path_y{m}(p),path_x{m}(p))==0) 
                power_rx(k,m) = power_rx(k,m) - loss_wall;
            end
        end
    end
    % Determine which Fog server (if any) is easiest to contact
    rxd_powers = reshape(power_rx(k,1:numfogsrv),numfogsrv,1);
    [max_power,im] = max(rxd_powers);
    if max_power < sensitivity_rx
        made_connection(k) = false;
        num_hops(k) = 0;
        next_hop(k) = 0;
    else
        made_connection(k) = true;
        num_hops(k) = 1;
        next_hop(k) = -im;
        if flag_plot
            plot(mote_col(k),mote_row(k),'bo');
            text(mote_col(k),mote_row(k),['  ',num2str(num_hops(k))]);
            plot(path_x{im},path_y{im},'c-'); 
        end
    end
end
num_connected = zeros(1,1);
num_connected(1) = nnz(made_connection);
% Continue as long as all motes haven't yet been connected and while 
% connections continue to be made
count_hops = 2;
while ~all(made_connection) && (num_connected(count_hops-1) > 0)
    num_connected(count_hops) = 0;
    for k = 1:nummotes
        if ~made_connection(k) % For each mote that hasn't yet been connected
            for m = 1:nummotes % Find closest mote that is already connected
                if (m==k) || ~made_connection(m) % If same mote or not connected
                    power_rx(k,m) = -inf;
                else % If connected, estimate distance between mote locations
                    diff_col(m) = abs(mote_col(k)-mote_col(m));
                    diff_row(m) = abs(mote_row(k)-mote_row(m));
                    dist_m(k,m) = pix2ft*ft2m*sqrt(diff_col(m)^2+diff_row(m)^2);
                    apathloss(k,m) = fspl(dist_m(k,m),c./freq);
                    power_rx(k,m) = power_tx - apathloss(k,m);
                    % Find straight path b/w motes, estimate #obstructions
                    npathpts(m) = max(diff_col(m),diff_row(m))+1;
                    fit_x = [mote_col(k),mote_col(m)]; 
                    fit_y = [mote_row(k),mote_row(m)];
                    pc = polyfit(fit_x,fit_y,1);
                    if (diff_col(m) >= diff_row(m))
                        path_x{m} = min(fit_x):max(fit_x);
                        path_y{m} = round(polyval(pc,path_x{m}));
                    else
                        if (min(fit_x) == max(fit_x))
                            path_x{m} = fit_x(1)+zeros(1,npathpts(m));
                            path_y{m} = min(fit_y):max(fit_y);
                        else
                        path_x{m} = linspace(min(fit_x),max(fit_x),npathpts(m));
                        path_y{m} = round(polyval(pc,path_x{m}));
                        path_x{m} = round(path_x{m});
                        end
                    end
                    for p = 1:npathpts(m)
                        if (I2(path_y{m}(p),path_x{m}(p))==0) 
                            power_rx(k,m) = power_rx(k,m) - loss_wall;
                        end
                    end
                end
            end
            rxd_powers = reshape(power_rx(k,1:nummotes),nummotes,1);
            [max_power(k),im(k)] = max(rxd_powers);
            if max_power(k) > sensitivity_rx
                make_connection(k) = true;
                num_hops(k) = count_hops;
                next_hop(k) = im(k);
                if flag_plot
                plot(mote_col(k),mote_row(k),'mo');
                text(mote_col(k),mote_row(k),['  ',num2str(num_hops(k))]);
                plot(path_x{im(k)},path_y{im(k)},'g-'); 
                end
                num_connected(count_hops) = num_connected(count_hops) + 1;
            end
        end
    end
    for k = 1:nummotes
        if make_connection(k)
            made_connection(k) = true;
            make_connection(k) = false;
        end
    end
    count_hops = count_hops + 1;
end
coverage = nnz(I3)/(nr*nc);
network_traffic = sum(num_hops);
percent_reachable = nnz(made_connection)/nummotes;
if (percent_reachable == 1)
    if flag_disp
    fprintf('Trial %d: Success! All %d sensor motes connected.',n,nummotes);
    end
else
    if flag_disp
    fprintf('Trial %d: Fail! Only %d of %d sensor motes connected.',n, ...
        nnz(made_connection),nummotes);
    end
    for k = 1:nummotes
        if ~made_connection(k)
            if flag_plot, plot(mote_col(k),mote_row(k),'rx'); end
        end
    end
end
if flag_plot
    figure(5); clf; imshow(I3); hold on; 
end
outputs(n).mote_row = mote_row;
outputs(n).mote_col = mote_col;
outputs(n).made_connection = made_connection;
outputs(n).num_hops = num_hops;
outputs(n).next_hop = next_hop;
outputs(n).coverage = coverage;
outputs(n).network_traffic = network_traffic;
outputs(n).percent_reachable = percent_reachable;
end % LOOP N: 1 to NUMTRIALS
save(sprintf('HICSS%02d.mat',nummotes),'outputs');
end % LOOP O: 1 to MAXNUMMOTES
end % FUNCTION HICSS3