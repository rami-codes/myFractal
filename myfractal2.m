% eegplugin_myfractal - Computes the fractal dimension of channels in an
% epoched EEG data set.
%
%
% Usage:
%   >> eegplugin_myfractal(EEG);                        % pop_up window
%   >> eegplugin_myfractal(EEG,method);                 %
%   >> eegplugin_myfractal(EEG,calculate_fd,method,channels,window size,window_overlap);


function [output_args] = eegplugin_myfractal2(EEG,calculate_fd,method,channels,window_size,window_overlap);

if nargin < 1
    help eegplugin_myfractal2
    error('You need to include atleast 1 argument when calling this function. Please read the help above.')
    return;
end
lastcom = [];
if nargin < 2
    popup = 1;
else
    popup = 0;
end
if nargin < 5
    window_size = 0;
    window_overlap = 0;
end

if popup
    [res userdata err structout] = inputgui( 'geometry', { [1 2] 1 [1 2] 1 [1 1] 1 [1 1] [1 1]}, ...
        'geomvert', [1 1 1 1 1 1 1 1], 'uilist', { ...
        {'style', 'text', 'string', 'Calculate FD?'}...
        { 'style', 'popupmenu', 'string', 'Yes|No' 'tag' 'calculate_fd' } { } {'style', 'text', 'string', 'Use FD Method:'}...
        { 'style', 'popupmenu', 'string', 'Calculate FD of each epoch then average|Average EEG epoch data then calculate FD|Concatenate all epochs then find FD|Window' 'tag' 'fd_method', 'userdata', 'continuous:off'}...
        { } {'style', 'text', 'string', 'Select channels to use ([] = all)'}   { 'Style', 'edit', 'string', '[]' 'tag' 'channels'} {}...
        {'style', 'text', 'string', 'Select size of window (powers of 2 only)'} { 'Style', 'edit', 'string', '64' 'tag' 'window_size'}...
        {'style', 'text', 'string', 'Select window overlap percent (0.0 to 1.0)'} { 'Style', 'edit', 'string', '0.50' 'tag' 'window_overlap', 'userdata', 'epoch:on'}}) ;
else
    structout.calculate_fd = calculate_fd;
    structout.method = method;
    structout.window_size = window_size;
    structout.fd_method = method;
    structout.channels = channels;
    structout.overlap = window_overlap;
end

folder_name = uigetdir;

if isempty(structout) == 0
    %if EEG.trials <=1
    %    error('You are using a continuous data set. You must use this extension with an EPOCHED dataset.')
    %end
    original_folder_name = pwd;
    if nargin < 2
        channels = eval( [ '[' structout.channels    ']' ] );
    end
    size_mat = size(EEG.data);
    calculate_fd = structout.calculate_fd;
    
    if isempty(channels)
        num_of_chans = EEG.nbchan;
        chans_list = [1:num_of_chans];
    else
        chans_list = channels;
        num_of_chans = length(chans_list);
    end
       
    if (calculate_fd == 0);
        if exist('fd_values') == 0
            disp('Please specify fd_values file directory');
            original_folder_name = pwd;
            load('fd_values');
        end
    end
end

window_increment = str2num(structout.window_overlap)*window_size;
if (str2num(structout.window_size) < window_increment)
    error('The size of the window must be greater than the window increment size. See documentation for more help.')
end

method_list = {'avg_fd' 'avg' 'concat' 'window'};
method = method_list{structout.fd_method};

if (length(size_mat) < 3)
    size_mat(3) = 1;
end

epochs = size_mat(3);

disp('Using the following channels: ')
disp(chans_list)
switch method
    case 'avg_fd'
        for i=1:size_mat(3) % epochs
            for k=1:num_of_chans % number of channels
                x.chan(chans_list(k)).epoch(i).fd = getFD(EEG.data(chans_list(k),:,i),'FDTool');
            end
        end
        
        %fd_values_mat = zeros(size_mat(3),size_mat(1));
        for i=1:size_mat(3) % epochs
            for k=1:num_of_chans % number of channels
                fd_values_mat(chans_list(k),i) = x.chan(chans_list(k)).epoch(i).fd ;
            end
        end
        
        %handles.A = figure(2) % create new figure
        
        %epoch_fds = size() compute size! finish me
        fd_values_mat(fd_values_mat==0) = [];
        
        % FIX ME
        %figure
        for k=1:size_mat(3) %epochs
            if (mod(k,37) == 0)
                %if (k>1)
                    figure
                %end
            end
            subplot(6,6,max(mod(k,37),1));
            topoplot(fd_values_mat(:,k)-mean(fd_values_mat(:,k)), EEG.chanlocs(chans_list));
            subplot_title = ['Epoch ' num2str(k)];
            title(subplot_title)
            epoch_fds(k,:,:) = fd_values_mat(:,k);
        end
        
        fd_values_avg = zeros(1,size_mat(1));
        for k=1:num_of_chans % channels
            fd_values_avg(chans_list(k)) = mean(fd_values_mat(chans_list(k),:));
        end
        
        fd_values_avg(fd_values_avg==0) = [];
        %handles.B = figure(3)
        figure
        topoplot(fd_values_avg-mean(fd_values_avg), EEG.chanlocs(chans_list));
        a = strcat('FD Averaged over Epochs: ',EEG.setname);
        title(a)
        colorbar()
        
        cd(folder_name)
        aa = strcat('FD values_',EEG.setname,'.mat');
        bb = strcat('FD',EEG.setname,'.fig');
        
        save(aa,'fd_values_avg')
        save('Epoch FD Values','epoch_fds')
        disp('Saving epoch data')
        saveas(gcf,bb);
        %save EPOCH DATA FIX ME
        
        cd(original_folder_name)

    case 'avg'
        for k=1:num_of_chans % channels
            y.chan(chans_list(k)).data = zeros(1,size_mat(2));
            for i=1:size_mat(3) % epochs
                for j=1:size_mat(2) %samples in epoch
                    y.chan(chans_list(k)).data(j) = y.chan(chans_list(k)).data(j) + EEG.data(chans_list(k),j,i);
                end
                
            end
        end
        for k=1:num_of_chans
            y.chan(chans_list(k)).data = y.chan(chans_list(k)).data/size_mat(3);
        end
        
        for k=1:num_of_chans
            y.chan(chans_list(k)).fd = getFD(y.chan(chans_list(k)).data,'FDTool');
            fd_values_avg(chans_list(k)) =  y.chan(chans_list(k)).fd;
        end
        
        fd_values_avg(fd_values_avg==0) = [];
        figure
        m=chans_list;
        bar(m,fd_values_avg-mean(fd_values_avg))
        a = strcat('FD of averaged data: ',EEG.setname);
        title(a)
        
        figure
        topo_data = fd_values_avg-mean(fd_values_avg);
        topoplot(topo_data, EEG.chanlocs(chans_list));
        title(a)
        colorbar()
        
        overall_mean = mean(fd_values_avg)
        use_subset = 0;
        if (use_subset == 1)
            figure
            bar(my_subset_array,fd_values_avg(my_subset_array)-mean(fd_values_avg(my_subset_array)))
            a = strcat('FD of Averaged EEG Epoch Data: ',EEG.setname);
            title(a)
            
            figure
            topo_data = fd_values_avg(my_subset_array)-mean(fd_values_avg(my_subset_array));
            topoplot(topo_data, EEG.chanlocs(my_subset_array));
            title(a)
        end
        
        cd(folder_name)
        disp('Select data save location: ');

        save('channel_fd_values','fd_values_avg');
        %save('fd_channel_means_averaged','fd_window_mean');
        
    case 'concat'
        for k=1:num_of_chans % number of channels
            x.chan(chans_list(k)).data = [];
            for i=1:size_mat(3) % number of epochs
                x.chan(chans_list(k)).data = [x.chan(chans_list(k)).data EEG.data(chans_list(k),:,i)];
            end
        end
        
        for k=1:num_of_chans % number of channels
            x.chan(chans_list(k)).fd = getFD(x.chan(chans_list(k)).data,'FDTool');
            fd_values_conc(chans_list(k)) =  x.chan(chans_list(k)).fd;
        end
        
        fd_values_conc(fd_values_conc==0) = [];
        
        m=1:num_of_chans;
        figure
        bar(m,fd_values_conc-mean(fd_values_conc))
        a = strcat('FD of concatenated data: ',EEG.setname);
        title(a)
        
        figure
        topo_data = fd_values_conc-mean(fd_values_conc);
        topoplot(topo_data, EEG.chanlocs(chans_list));
        title(a)
        colorbar()
        
        cd(folder_name)
        %a2 = strcat(a,'.mat')
        save('FD_concatenated_data.mat','fd_values_conc');
        
    case 'window'
        disp('This method is exceptionally computationally intensive! Please be patient.')
     
        window_size = 2^(nextpow2(str2num(structout.window_size))); % find next power of 2 given window_size
        window_increment = str2num(structout.window_overlap)*window_size;
        
        back(1) = 1;
        front(1) = window_size;
        
        num_windows = 0;
        for k=2:(size_mat(2)) 
            if ((front(k-1) + window_increment) < size_mat(2)) % ensures front is within data range
                back(k) = back(k-1) + window_increment;
                front(k) = front(k-1) + window_increment;
                num_windows = num_windows + 1;
            else 
                break
            end
        end
        
        for n=1:size_mat(3) % number of epochs
            for j=1:num_of_chans % number of channels
                for i=1:(num_windows) % number of windows
                    fd_values.chan(chans_list(j)).epoch(n).window(i) = getFD(EEG.data(chans_list(j),back(i):front(i),n),'FDTool');
                    fd_values_mat(chans_list(j),i,n) = fd_values.chan(chans_list(j)).epoch(n).window(i);
                    fd_mean_mat(chans_list(j),i) = mean(fd_values_mat(chans_list(j),i,:));
                end
            end
        end
        
       
        for j=1:num_of_chans %num of channels
            fd_window_mean(j) = mean(fd_mean_mat(j,:));
        end
       
        %h=1:length(fd_mean_mat);
        h = linspace(0,length(EEG.times)/EEG.srate,length(fd_mean_mat));
        figure
        plot(h,fd_mean_mat')
        xlim([0 length(EEG.times)/EEG.srate])
        a = strcat('Windowed FD of Epoched Data: ',EEG.setname);
        title(a)
        xlabel('Time')
        ylabel('Fractal Dimension')
        
        
        figure
        topoplot(fd_window_mean-mean(fd_window_mean), EEG.chanlocs(chans_list));
        title('FD Topoplot')
        colormap('default')
        colorbar()
        
        disp('Saving data matrices in folder');
        save('fd_channel_means_overtime.mat','fd_mean_mat');
        save('fd_channel_means_averaged.mat','fd_window_mean');
        
        
end
disp('Done!');
end
