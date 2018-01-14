function [ cs, cb, cr ] = pres_ui( plot_frame, n, dur, tps, maxfps )
%PRES_UI Generate user interface on current figure.
%   Args:
%     plot_frame    Handle of the plot-updating function. The passed
%                     arguments are `i` and `frames` such that `frames(i)`
%                     is the current frame number of the animation
%                     with respect to 1:n.
%     n             Total number of frames.
%     dur           Animation duration in simulation time units.
%     tps           Animation rate -- simulation time units per second.
%     maxfps        Framerate limit (default: 30).
%
%   Returns:
%     cs, cb, cr    Widget controls -- slider, animate button and reset
%                   button respectively.

nin = nargin;
if nin < 5
    [fps, frames] = anim_tps2fps(tps, dur, n);
else
    [fps, frames] = anim_tps2fps(tps, dur, n, maxfps);
end
L = length(frames);
sstep = 0.04;

fig = gcf;
cs = uicontrol(fig, 'Style', 'Slider',...
    'String', 'Time',...
    'Min', 1, 'Max', L, 'Value', 1,...
    'SliderStep', [1/L, max([sstep,1/L])],...
    'Position', [20, 50, 200, 20]);
ct = uicontrol(fig, 'Style', 'Edit',...
    'String', tps, 'HorizontalAlignment', 'right',...
    'Position', [160, 20, 60, 20]);
cb = uicontrol(fig, 'Style', 'pushbutton',...
    'String', 'Animate',...
    'Position', [20,20,100,20]);
cr = uicontrol(fig, 'Style', 'pushbutton',...
    'String', 'Reset',...
    'Position', [120,20,60,20]);

state = false;
i = 1;

    function ui_state_on()
        state = true;
        cb.String = 'Pause';
    end

    function ui_state_off()
        state = false;
        cb.String = 'Animate';
    end

    function ui_slide(src,~)
        ui_state_off;
        i = round(src.Value);
        plot_frame(i, frames);
    end

    function ui_anim(~, ~)
        if ~state
            if i == L
                i = 1;
            end
            ui_state_on;
            while i <= L && state
                cs.Value = i;
                plot_frame(i,frames)
                pause(1/fps);
                i = i+1;
            end
            if i == L+1
                ui_state_off;
                i = L;
            end
        else
            ui_state_off;
        end
    end

    function ui_reset(~,~)
        ui_state_off
        i = 1;
        plot_frame(i,frames);
        cs.Value = 1;
    end

    function ui_change_tps(src,~)
        cfr = frames(i);
        new_tps = str2double(src.String);
        if isnan(new_tps)
            warning('textbox: tps is not a number!');
        else
            if nin < 5
                [fps, frames] = anim_tps2fps(new_tps, dur, n);
            else
                [fps, frames] = anim_tps2fps(new_tps, dur, n, maxfps);
            end
            L = length(frames);
            if i <= L
                i = find(frames >= cfr, 1);
                cs.Max = L;
                cs.Value = i;
            else
                i = find(frames >= cfr, 1);
                cs.Value = i;
                cs.Max = L;
            end
            cs.SliderStep = [1/L, max([sstep,1/L])];
        end
    end

cs.Callback = @ui_slide;
ct.Callback = @ui_change_tps;
cb.Callback = @ui_anim;
cr.Callback = @ui_reset;

end

