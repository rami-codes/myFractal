% eegplugin_myfractal - Computes the fractal dimension of channels in an
% epoched EEG data set.
%
%
% Usage:
%   >> eegplugin_myfractal(EEG);                        % pop_up window
%   >> eegplugin_myfractal(EEG,method);                 %
%   >> eegplugin_myfractal(EEG,method,window size);


function vers = eegplugin_myfractal2(fig,try_strings,catch_strings)

vers = 'my_fractal 2.0';

menu = findobj(fig,'Tag','tools');
%abovemenupos = get(findobj(toolsmenu,'Label','Remove components'),'position');


uimenu( menu, 'label', 'myFractal (Compute FD)', 'callback',...
    [ 'myfractal2(EEG);' ]); 

%,'position',abovemenupos+1