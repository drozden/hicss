function [I2,nr,nc] = fix_gmap(filename,fogsrv_col,fogsrv_row,flag_plot)
if nargin < 4
    flag_plot = true;
    if nargin < 3
        fogsrv_row = [236,330];
        if nargin < 2
            fogsrv_col = [1031,1324];
            if nargin < 1
                filename = 'area3.JPG';
            end
        end
    end
end
I0 = imread(filename);
figure(1); clf; imshow(I0); hold on; 
for ii = 1:length(fogsrv_col)
    plot(fogsrv_col(ii),fogsrv_row(ii),'mp');
end
[nr,nc,~] = size(I0);
I = double(I0);
I1 = rgb2gray(I0); 
I2 = zeros(size(I1),class(I1))+128;
% Colors in Bogard Hall
col_bgrd = [241,242,236; 242,239,232; 242,241,234; 242,241,237; 242,242,234; ...
    243,240,233; 243,240,235; 244,240,237; 244,241,234; ...
    245,242,237; 245,244,240; 246,243,238; 249,244,238];
% Identify locations of all ground/road/lot pixels 
nb = 3;
for ii = 1:nr
    for jj = 1:nc
        this_col = reshape(I(ii,jj,:),1,3);
        if all(this_col>=250) % Road
            I(ii,jj,1:3) = [255,255,255];
            I2(max(1,ii-nb):min(nr,ii+nb),max(1,jj-nb):min(nc,jj+nb)) = 255;
        elseif isequal(this_col,[242,243,245]) % Ground
            I2(max(1,ii-nb):min(nr,ii+nb),max(1,jj-nb):min(nc,jj+nb)) = 255;
        elseif isequal(this_col,[243,243,243]) % Ground
            I2(max(1,ii-nb):min(nr,ii+nb),max(1,jj-nb):min(nc,jj+nb)) = 255;
        elseif isequal(this_col,[237,241,242]) % Ground
            I2(max(1,ii-nb):min(nr,ii+nb),max(1,jj-nb):min(nc,jj+nb)) = 255;
        elseif isequal(this_col,[239,243,242]) % Ground
            I2(max(1,ii-nb):min(nr,ii+nb),max(1,jj-nb):min(nc,jj+nb)) = 255;
        elseif isequal(this_col,[248,249,251]) % Parking Lot
            I2(max(1,ii-nb):min(nr,ii+nb),max(1,jj-nb):min(nc,jj+nb)) = 255;
        elseif isequal(this_col,[251,246,240]) % Bldg: Book
            I2(max(1,ii-nb):min(nr,ii+nb),max(1,jj-nb):min(nc,jj+nb)) = 0;
        elseif isequal(this_col,[251,246,242]) % Bldg: PO
            I2(max(1,ii-nb):min(nr,ii+nb),max(1,jj-nb):min(nc,jj+nb)) = 0;
        elseif isequal(this_col,[255,246,237]) % Bldg: Cafe
            I2(max(1,ii-nb):min(nr,ii+nb),max(1,jj-nb):min(nc,jj+nb)) = 0;
        end
    end
end
ns = 3; thr255 = 6;
if flag_plot
    figure(2); clf; imshow(I2); hold on; 
    for ii = 1:length(fogsrv_col)
        plot(fogsrv_col(ii),fogsrv_row(ii),'mp');
    end
end
for ii = 1:nr
    for jj = 1:nc
        if (I2(ii,jj) ~= 255)
        this_col = reshape(I(ii,jj,:),1,3);
        if isequal(this_col,[237,237,237]) % Bldg: Neth
            nh = I2(max(1,ii-ns):min(nr,ii+ns),max(1,jj-ns):min(nc,jj+ns));
            n255 = nnz(nh==255);
            if (n255 < thr255)
                I2(max(1,ii-nb):min(nr,ii+nb),max(1,jj-nb):min(nc,jj+nb))=0;
            end
        elseif isequal(this_col,[239,239,239]) % Bldg: Libr
            nh = I2(max(1,ii-ns):min(nr,ii+ns),max(1,jj-ns):min(nc,jj+ns));
            n255 = nnz(nh==255);
            if (n255 < thr255)
                I2(max(1,ii-nb):min(nr,ii+nb),max(1,jj-nb):min(nc,jj+nb))=0;
            end
        elseif ismember(this_col, col_bgrd) % Bldg: Bogard
            nh = I2(max(1,ii-ns):min(nr,ii+ns),max(1,jj-ns):min(nc,jj+ns));
            n255 = nnz(nh==255);
            if (n255 < thr255), I2(ii,jj) = 0; end
        end
        end
    end
end
if flag_plot
    figure(3); clf; imshow(I2); hold on; 
    for ii = 1:length(fogsrv_col)
        plot(fogsrv_col(ii),fogsrv_row(ii),'mp');
    end
end
% For remaining, unidentified pixels, find most likely pixel
% using pixels in the neighborhood
ns = 3; cenh = 0;
for ii = 1:nr
    for jj = 1:nc
        this_I = I2(ii,jj);
        if (this_I == 128)
            nh = I2(max(1,ii-ns):min(nr,ii+ns),max(1,jj-ns):min(nc,jj+ns));
            n0 = nnz(nh==0);
            n255 = nnz(nh==255);
            if (n0>n255)
                I2(ii,jj) = 0;
            elseif (n255>=n0)
                I2(ii,jj) = 255;
            else
                cenh = cenh + 1;
                plot(jj,ii,'rx');
                disp(['Equal neighborhood: ',num2str(cenh)]);
            end
        end
    end
end
if flag_plot
    figure(4); clf; imshow(I2); hold on; 
    for ii = 1:length(fogsrv_col)
        plot(fogsrv_col(ii),fogsrv_row(ii),'mp');
    end
end
end % FUNCTION FIX_GMAP