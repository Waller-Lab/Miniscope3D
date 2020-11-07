% Read in PSF
% This is designed to work with measurements in a folder, then in a
% parallel folder, save the recons (i.e. measurements../recons/)



%psf_path = 'D:\Kyrollos\RandoscopeNanoscribe\RandoscopeNanoscribe\Miniscope3D\psf_svd_12comps_23z_240xy_20190619';
%psf_path = 'D:\Antipa\Randoscopev2_PSFs\Data_8_21_2019\SVD_2_5um_PSF_20um_1';

psf_path = 'D:\Antipa\Randoscopev2_PSFs\20190912_recalibration\SVD_2p5_um_PSF_5um_1_green_channel';
%psf_path = 'T:\Antipa\Randoscopev2_PSFs\20190912_recalibration\SVD_2p5_um_PSF_5um_1_green_channel';
comps_path = [psf_path,'\SVD_2_5um_PSF_5um_1_ds2_components_green_SubAvg.mat'];
weights_path = [psf_path,'\SVD_2_5um_PSF_5um_1_ds2_weights_interp_green_SubAvg.mat'];

%%
%comps_path = [psf_path,'/SVD_2_5um_PSF_5um_1_ds2_components_green_NoFro.mat'];
%weights_path = [psf_path,'/SVD_2_5um_PSF_5um_1_ds2_weights_interp_green_NoFro.mat'];
fprintf('loading components\n')

h_in = load(comps_path);
fprintf('done.\nLoading weights\n')
weights_in = load(weights_path);
fprintf('done loading PSF data\n')
%%

%Get names of files/paths
[meas_name,data_path,~] = uigetfile('*.*','Select measurement','T:\Randoscope\RandoscopeV2_data');
dots = strfind(meas_name,'.');
fext = meas_name(dots(end):end);
if strcmpi(fext,'.tif')
    [bg_name, bg_path,~] = uigetfile('*.*',['Select background for ',meas_name],fullfile([data_path,'../']));
else
    bg_name = 'NONE';
    bg_path = 'NONE';
end

meas_path = [data_path,meas_name];
bg_path = [bg_path,bg_name];

if strcmpi(fext,'.tif')
    ome = strfind(meas_name,'.ome');
    bg_ome = strfind(bg_name,'.ome');
    meta_name = [meas_name(1:ome-1),'_metadata.txt'];
    bg_meta_name = [bg_name(1:bg_ome-1),'_metadata.txt'];

    % Inline function to open a metadata file, convert to characters then
    % decode (it's in json format).
    parse_json = @(x)jsondecode(transpose(fread(fopen(x),'*char')));

    % Construct full system paths to images
    meta_path = [data_path,meta_name];
    bg_meta_path = [bg_path,bg_meta_name];


    % Get json data from measurement and background images

    file_info = parse_json(meta_path);
    bg_info = parse_json(bg_meta_path);

    params.data_format = file_info.FrameKey_0_0_0.x50890959_DataFormat;
    params.bg_format = bg_info.FrameKey_0_0_0.x50890959_DataFormat;
    params.ds_raw = file_info.FrameKey_0_0_0.Binning;
    params.ds_bg = bg_info.FrameKey_0_0_0.Binning;
elseif strcmpi(fext,'.mat')
    params.data_format = 'mat';
    params.bg_format = 'NONE';
    params.ds_raw = 4;
    mat_var_name = 'vid_bgrm_ds';
    meta_path = 'NONE';
    file_info='NONE';
    bg_info='NONE';
end

        


params.demosaic = contains(lower(params.data_format),'raw');
params.demosaic_bg = contains(lower(params.bg_format),'raw');
params.bg_name = bg_name;
params.meas_path = meas_path;
params.meta_path = meta_path;
params.bg_path = bg_path;
params.meas_info = file_info;
params.bg_info = bg_info;
%%
%Waterbear_20190905\waterbear_big_lastone_20_3_30ms';
%bg_path = 'D:\Randoscope\RandoscopeV2_data\Waterbear_20190905\waterbear_big_lastone_bck_20_3_30ms_1';
 %= 'waterbear_big_lastone_3_MMStack_Default.ome.tif';
%bg_name = 'waterbear_big_lastone_bck_20_3_30ms_1_MMStack_Default.ome.tif';



    

%for zd = 9
    
    %data_path = 'Z:\kyrollos\RandoscopeNanoscribe\Nanoscribe_pdms\Data_8_21_2019\real_res_target_10um_1';   %<--folder where the measurements are
    %bg_path = 'Z:\kyrollos\RandoscopeNanoscribe\Nanoscribe_pdms\Data_8_21_2019\bck_real_res_target_10um_1';
   
    
    params.data_tiff_format = 'time';   %Use 'time' if tiff stacks are at the same location over time, use 'z' if they are z stacks'
    params.tiff_color = 2;    %use 'rgb' or 'mono'. Use number (1,2,3) for r,g, or b only
    params.meas_depth = 82;    %If using 3D tiff or list of files, which slice was processed?
    
    params.ds_z = 1;   %z downsampling ratio
    params.meas_bias = 0;
   
    params.ds = 4;  % Global downsampling ratio (i.e.final image-to-sensor ratio)
    params.ds_psf = 2;   %PSf downsample ratio (how much to further downsample -- if preprocessing included downsampling, use 1)
    params.ds_meas = params.ds/params.ds_raw;   % How much to further downsample measurement?
    params.z_range = 1:44; %Must be even number!! Range of z slices to be solved for. If this is a scalar, 2D. Use this for subsampling z also (e.g. 1:4:... to do every 4th image)
    params.rank = 12;
    useGpu = 1; %cannot fit ds=2 on gpu unless we limit z range!!!!
    params.psf_norm = 'fro';   %Use max, slice, fro, or none
    
    %meas_name = ['real_res_target_10um_1_MMStack_Img_',num2str(params.meas_depth),'_000_000.ome.tif'];    %<--- name of measurement
    %bg_name = ['bck_real_res_target_10um_1_MMStack_Img_',num2str(params.meas_depth),'_000_000.ome.tif'];
  
    

    
    
    % Make sure h and weights are in order y,x,z,rank
    fprintf('permuting PSF data\n')
    h = permute(h_in.comps_out(:,:,1:params.rank,params.z_range),[1,2,4,3]);
    weights = permute(weights_in.weights_out(:,:,1:params.rank,params.z_range),[1,2,4,3]);
    fprintf('Done permuting. Resampling PSF\n');
    
    %clear h_in;
    %clear weights_in;
    h = single(imresize(squeeze(h),1/params.ds_psf,'box'));
    weights = single(imresize(squeeze(weights),1/params.ds_psf,'box'));
    
    % Normalize weights to have maximum sum through rank of 1
    weights_norm = max(sum(weights(size(weights,1)/2,size(weights,2)/2,:,:),4),[],3);
    weights = weights/weights_norm;
    fprintf('Done. PSF ready!\n')
    %clear h_permute;
    %clear weights_permute;

    %%
    switch lower(params.psf_norm)
        case('max')
            h = h/max(h(:));
        case('none')
        case('fro')
            h = h/norm(vec(h));
        case('slice')
            for sl = 1:Nz
                slice_norm = norm(h(:,:,sl,1),'fro');
                for cp = 1:Nr
                    h(:,:,sl,cp) = h(:,:,sl,cp)/slice_norm;
                end
            end
    end
    
    H = fft2(ifftshift(ifftshift(h,1),2));
    Hconj = conj(H);
    if useGpu
        H = gpuArray(H);
        Hconj = gpuArray(Hconj);
        weights = gpuArray(weights);
    end
    
    
    % Read in data
    
%%
    init_style = 'loaded';   %Use 'loaded' to load initialization, 'zeros' to start from scratch. Admm will run 2D deconv, then replicate result to all time points
    im_tag = 'rank1_soft_of_TV_TEST';
    if strcmpi(params.data_format,'mat')
        data_in = load(meas_path,mat_var_name);
        data_in = data_in.(mat_var_name);
    end
        
        
    for meas_slice = 1:33
        
        params.meas_slice = meas_slice;   %Slices to load from tiff stack. If 'all' used, it will average.
        %params.meas_slice = 'all';
        if ~strcmpi(params.data_format,'mat')
            switch lower(params.data_tiff_format)
                case('z')
                    data_raw = double(read_tiff_stack(meas_path,params.ds_meas,params.meas_depth));
                    bg_in =  double(read_tiff_stack(bg_path,params.ds_meas,params.meas_depth));
                case('time')
                    bg_raw = read_tiff_stack(bg_path,1);
                    if strcmpi(params.meas_slice,'all')
                        data_raw = mean(double(read_tiff_stack(meas_path,1)),4);   %Average out the time variable
                        if params.demosaic
                            data_demos = imresize(double(demosaic(uint16(data_in),'grbg')),params.ds_raw/params.ds,'box');
                        else
                            data_demos = imresize(data_raw,params.ds_raw/params.ds,'box');
                        end
                    else
                        data_in = read_tiff_stack(meas_path,1,params.meas_slice);


                        if params.demosaic
                            data_demos = imresize(double(demosaic(uint16(data_in),'grbg')),params.ds_raw/params.ds,'box');
                            %grbg
                        else
                            data_demos = imresize(mean(data_in,4),params.ds_raw/params.ds,'box');
                        end
                    end
                    if params.demosaic_bg
                        bg_in = imresize(double(demosaic(uint16(mean(bg_raw,3)),'grbg')),params.ds_bg/params.ds,'box');
                    else
                        bg_in = imresize(mean(double(bg_raw),4),params.ds_bg/params.ds,'box');
                    end

                    % data_raw = data_raw(:,:,:,1);

            end
            
            if strcmpi(params.tiff_color,'rgb')
                data = mean(data_demos,3);
                bg = mean(bg_in,3);   %Average out color. Change to (:,:,color) to select one channel
            elseif isnumeric(params.tiff_color)
                data = data_demos(:,:,params.tiff_color);
                bg = bg_in(:,:,params.tiff_color);
            end
        else
            data = data_in(:,:,meas_slice);
            bg = 0;
        end
        
      
        data = data - bg - params.meas_bias;
        b = data/max(data(:));
        
        
        % data_r = data_in(:,:,1);
        % data_g = data_in(:,:,2);
        % data_b = data_in(:,:,3);
        
        
        %Nx = size(h,2);
        %Ny = size(h,1);
        if numel(size(h)) == 3
            [Ny, Nx, Nr] = size(h);
            Nz = 1;
        else
            [Ny, Nx, Nz, Nr] = size(h);
        end
        
        %define crop and pad operators to handle 2D fft convolution
        pad2d = @(x)padarray(x,[size(h,1)/2,size(h,2)/2],0,'both');
        ccL = size(h,2)/2+1;
        ccU = 3*size(h,2)/2;
        rcL = size(h,1)/2+1;
        rcU = 3*size(h,1)/2;
        
        %cc = gpuArray((size(h,2)/2+1):(3*size(h,2)/2));
        %rc = gpuArray((size(h,1)/2+1):(3*size(h,1)/2));
        crop2d = @(x)x(rcL:rcU,ccL:ccU);
        
        
        
        
        
        
        
        if strcmpi(init_style, 'zeros')
            xinit = zeros(Ny, Nx, Nz);
        elseif strcmpi(init_style,'loaded')
            if ~exist('xinit')
                xinit = zeros(Ny,Nx,Nz);
            else
                xinit = xhat_out(:,:,:);
            end
        elseif strcmpi(init_style,'admm')
            xinit_2d = gpuArray(single(zeros(Ny, Nx, 3)));
            
            for n = 1:3
                xinit_2d(:,:,n) = admm2d_solver(gpuArray(single(b(:,:,n))), gpuArray(single(h(:,:,n))),[],.001);
                
                imagesc(2*xinit_2d/max(xinit_2d(:)))
            end
        end
        
        
        
        
        
        
        
        options.color_map = 'parula';
        
        
        
        options.convTol = 15e-4;
        
        %options.xsize = [256,256];
        options.maxIter = 2000;
        options.residTol = 5e-5;
        options.momentum = 'nesterov';
        options.disp_figs = 1;
        options.disp_fig_interval = 40;   %display image this often
        if Nz == 1
            options.xsize = [Ny, Nx];
        else
            options.xsize=[Ny, Nx, Nz];
        end
        options.print_interval = 20;
        
        figure(2)
        clf
        imagesc(b)
        axis image
        
        h1 = figure(1);
        clf
        options.fighandle = h1;
        nocrop = @(x)x;
        options.known_input = 0;
        
        
        
        

        large = 0;
        if Nz > 1
            if large == 0
                A = @(x)A_svd_3d(x, weights,H);

                Aadj = @(y)A_adj_svd_3d(y, weights, Hconj);
            else
                weights=gpuArray(weights);
                H = gpuArray(H);
                Hconj = gpuArray(Hconj);
                b = gpuArray(single(b));
                A = @(x)A_svd_3d_large(x,weights,H);
                Aadj = @(y)A_adj_svd_3d_large(y, weights, Hconj);
            end
        elseif Nz == 1
            A = @(x)A_svd(H, weights, x, nocrop);
            Aadj = @(y)A_adj_svd(Hconj,weights,y,nocrop);
        end
        
        
        
        
               %options.stepsize = .1e-3; for ds=4
        if params.ds == 4
            if strcmpi(params.psf_norm ,'fro')
                if Nz == 18
                    options.stepsize = 3e-3;
                elseif Nz == 12
                    options.stepsize = .4e-2;
                    fprintf('foo\n')
                elseif Nz == 14
                    options.stepsize = 4e-3;
                elseif Nz == 20
                    if params.rank  == 12
                        options.stepsize = 3e-3;
                    elseif params.rank == 8
                        options.stepsize = 1e-3;
                    elseif params.rank == 18
                        options.stepsize = 4e-3;
                    
                       
                    end
                elseif Nz>20
                    options.stepsize = .014;  %015 is nice?
                end
                
            else
                options.stepsize = 3e-6;
            end
            
        elseif params.ds == 2
            options.stepsize = 0.7e-3;
        end
        %.3e-4 is good for waterbears
        params.tau1 = options.stepsize*.3e-3; %was 0.5e-7   %.000005 works pretty well for v1 camera, .0002 for v2
        params.tau_soft = options.stepsize * 3e-3;
        tau_iso = (.25e-4);
        params.z_tv_weight = 1;    %z weighting in anisotropic TV
        tau2 = .001  %Auxilliary
        TVnorm3d = @(x)sum(sum(sum(abs(x))));
        
        
        if useGpu
      
            grad_handle = @(x)linear_gradient_b(x, A, Aadj, gpuArray(single(b)));
   
            params.tau1 = gpuArray(params.tau1);
            params.tau_soft = gpuArray(params.tau_soft);
            tau_iso = gpuArray(tau_iso);
            params.z_tv_weight = gpuArray(params.z_tv_weight);
            options.stepsize = gpuArray(options.stepsize);
            
        else
            if ~large
                grad_handle = @(x)linear_gradient_b(x, A, Aadj, single(b));
            else
                grad_handle = @(x)linear_gradient_large(x,A,Aadj,gpuArray(single(b)));
            end
            
        end
        
        %Prox
        %prox_handle = @(x)deal(x.*(x>=0), abs(sum(sum(sum(x(x<0))))));
        
        
       
        %prox_handle = @(x)deal(1/3*(x.*(x>=0) + soft(x, tau2) + tv3dApproxHaar(x, params.tau1)), TVnorm3d(x));
        
        
 
        
        if ~strcmpi(params.data_format,'mat')
            if Nz>1
                prox_handle = @(x)deal(1/2*(max(x,0) + (tv3d_iso_Haar((x), params.tau1, params.z_tv_weight))), params.tau1*TVnorm3d(x));
            elseif Nz == 1
                prox_handle = @(x)deal(.5*tv2d_aniso_haar(x,params.tau1*options.stepsize) + ...
                    .5*max(x,0), params.tau1*options.stepsize*TVnorm(x));
            end
        else
%             prox_handle = @(x)deal(.5*(soft(x,params.tau_soft) +...
%                 tv3d_iso_Haar(x, params.tau1, params.z_tv_weight)), ...
%                 params.tau1*TVnorm3d(x));
            prox_handle = @(x)deal(soft(tv3d_iso_Haar(x, params.tau1, params.z_tv_weight),params.tau_soft), ...
                params.tau1*TVnorm3d(x));
            %prox_handle=@(x)deal(soft(x,params.tau_soft),params.tau_soft*sum(abs(vec(x))));
        end
        TVpars.epsilon = 1e-7;
        TVpars.MAXITER = 100;
        TVpars.alpha = .3;
        %prox_handle = @(x)deal(hsvid_TV3DFista(x, tau_iso, 0, 10, TVpars) , hsvid_TVnorm3d(x));
        
        if strcmpi(init_style, 'zeros')
            xinit = zeros(Ny, Nx, Nz);
            
        end
        
        if useGpu
            
            TVpars.epsilon = gpuArray(TVpars.epsilon);
            TVpars.MAXITER = gpuArray(TVpars.MAXITER);
            TVpars.alpha = gpuArray(TVpars.alpha);
            xinit = gpuArray(single(xinit));
            success = false;
            while success == false   %This shouldn't be necessary, but it deals with restarting when GPU runs OOM
                try
                    [xhat, f2] = MiniscopeFISTA(grad_handle,prox_handle,xinit,options);
                    success = true;
                catch
                    success = false;
                end
            end
                    
        else
            if large
                xinit = gpuArray(xinit);
            end
            [xhat, f2] = MiniscopeFISTA(grad_handle,prox_handle,xinit,options);
        end
        
        
        
        
        
        
        
        datestamp = datetime;
        tiff_string = sprintf('%03d',meas_slice);
        date_string = datestr(datestamp,'yyyy-mmm-dd_HHMMSS');
        save_str = ['../recons/',date_string,'_',meas_name(1:end-4),'_',im_tag,'_',tiff_string];
        full_path = fullfile(data_path,save_str);
        mkdir(full_path);
        
        
        
        imout = gather(xhat/prctile(xhat(:),100*(numel(xhat)-10)/numel(xhat)));   %Saturate only 10 pixels
        xhat_out = gather(xhat);
        
        params.tau1 = gather(params.tau1);
        params.tau_soft = gather(params.tau_soft);
        imbase = meas_name(1:end-4);
        mkdir([full_path, '/png/']);
        filebase = [full_path, '/png/', imbase];
        f_out = gather(f2);
        out_names = {};
        for n= 1:size(imout,3)
            out_names{n} = [filebase,'_',sprintf('Z_%.3i_T_',params.z_range(n)),...
                tiff_string,'_',im_tag,'.png'];
            imwrite(imout(:,:,n),out_names{n});
            fprintf('writing image %i of %i\n',n,size(xhat,3))
        end
        
        fprintf('zipping...\n')
        zip([full_path, '/png/', imbase],out_names)
        fprintf('done zipping\n')
        
        fprintf('writing .mat\n')
        options.fighandle = []
        
        save([full_path,'/',meas_name(1:end-4),'_',date_string,'_',im_tag,'_',tiff_string,'.mat'], 'tau_iso','TVpars','xhat_out', 'options', 'comps_path','weights_path', 'b','params')
        fprintf('done writing .mat\n')
        % gpuDevice(1)
        clear xhat
        clear f2


        if params.ds == 2
            gpuDevice(1)
        end
    end
%end

%%
% imagesc(brain_recon.xhat(:,:,10))
% axis image