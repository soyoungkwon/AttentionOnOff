% function DotDemo
% dot motion demo using SCREEN('DrawDots') subfunction
% Basic : Full screen dot motion.
% 1. Dot move in certain direction (Exp & CW/Red &CW/Exp & CCW/Red & CCW) (with noise, moves random dir)
%  & sometimes the direction change (e.g., EXp & CW --> Red & CW)
% 2. Color can be either R/G/B (with noise)
%  & sometimes the color change (e.g., Red --> Green)
%
% Prerequisite: Psychtoolbox 3.0.8
% fname_cols =  'colours_lin.mat';
% fname_LUT = 'spc_mrz3t_gammaLUT.txt';
% img_gammaConvert, wait_MRtrigger_multi.m
% original code /Applications/Psychtoolbox/PsychDemos/DotDemo.m by Keith Schneider, 2004/12/13
% modified code:
%  2009/12/08 direction U/D/L/R  by Soyoung Kwon
%  2010/01/20 direction Clockwise
clear all;
eyetrack =1;
auto = 1;
if auto == 1;
    stimuli_pc = input('stimuli pc?? yes(1) or no(0) : (due to screen number)');
    scan= input('scanning?? scan(1) or not scan(0) :');
    subject=input('subject name: ','s');
    sess=(input('which task? (psych: psychophysics, short, sess1: C/M/R, sess3: R/M/C, sess4: M/C/R)  ','s'));
    os = input('whic os? win, mac ? ', 's');
else
    stimuli_pc = 1;
    scan = 1;
    subject = 'test';% '19berhan' ;%'michael';%'marie';
    load runindex runindex
    sess = ['sess', num2str(runindex)];
    os = 'mac';%'win' ;%mac';%'win';
end
refreshrate= 60;
short_atten = 20; %in sec (20)
long_atten = 120; % in sec (120)
sequ = repmat([1 0], 1,8 ); %% total 16 times!


if strcmpi('psych', sess)
    attention_frame = short_atten*refreshrate;
    sequ = repmat([1 0], 1, 8);
    n = size(sequ,2);
    p_dir = 0.9;
elseif strcmpi('short', sess)
    attention_frame = short_atten*refreshrate; % each attention time duration 60 --> refresh rate
    sequ = repmat([1 0], 1,8 );
    n = size(sequ,2);
elseif strmatch('sess', sess) % localizer, 20s
    attention_frame = long_atten*refreshrate; % each attention time duration 60 --> refresh rate
    if strcmpi('sess1', sess) || strcmpi('sess3', sess) || strcmpi('sess5', sess) || strcmpi('sess7', sess) 
        sequ = repmat([1 0], 1,2);
    elseif strcmpi('sess2', sess) || strcmpi('sess4', sess) || strcmpi('sess6', sess) || strcmpi('sess8', sess) 
        sequ = repmat([0 1], 1, 2);
    end
    n = size(sequ,2);
end

%% p_dir!!
if strcmpi('psych', sess)
else
    p_dir = 0.65;
end
close all; clc;



%% -----------set dot field parameters-------------

% nframes     = 60*60*3*4; % number of animation frames in loop
mon_width   = 39;   % horizontal dimension of viewable screen (cm)
v_dist      = 60;   % viewing distance (cm)
max_d       = 15;   % maximum radius of  annulus (degrees)
min_d       = 1;    % minumum
dot_w       = 0.1;  % width of dot (deg)
fix_r       = 0.15; % radius of fixation point (deg)
f_kill      = 0.05; % fraction of dots to kill each frame (limited lifetime)
differentcolors =1; % Use a different color for each point if == 1. Use common color white if == 0.
differentsizes = 2; % Use different sizes for each point if >= 1. Use one common size if == 0.
waitframes = 1;     % Show new dot-images at each waitframes'th monitor refresh.

% Quest
% p =1;
total_dots = 600;
ndots = round(total_dots * 1);
%% difficulty!!!!!!
% how often changes! how much later detect counts as resp=1
color_change_rate = 60;
color_change_variance =10;
color_change_det_th = 60;
color_change_det_min = 10;
dir_change_rate = 60;
dir_change_variance = 10;
dir_change_det_th = 80;
dir_change_det_min = 20;
% proportion of one dir/col dots
if strcmpi(sess, 'psych')
    p_dir = 0.8;
%     p_col = 0.8;
    tGuessSd= 0.3;
else
    p_dir=floor(p_dir*100)/100;%0.6;%0.8;
    tGuessSd=0.05;
end
lum = 255*2;
% p_col_array =[]; 
p_dir_array = [];
tGuess_dir = log10(p_dir);
% tGuessSd = 0.05; %0.3-initial
pThreshold = 0.8; beta = 3.5; delta = 0.01; gamma = 0.05; grain = 0.05;
% q_col = QuestCreate(tGuess_col,tGuessSd,pThreshold,beta,delta,gamma,grain);
q_dir = QuestCreate(tGuess_dir,tGuessSd,pThreshold,beta,delta,gamma,grain);

% difficulty control
dot_speed   = 3;    % dot speed (deg/sec)
% dt = 0.01;    % theta change constant
% dt_tm = 0.01;

%% --------- open the screen--------------------------
AssertOpenGL;
doublebuffer=1
screens=Screen('Screens');
if stimuli_pc,
    screenNumber = 1;
else
    screenNumber=max(screens);
end

%% ------- keyboard setup ---------------
KbName('UnifyKeyNames');
escapeKey = KbName('Escape');
leftKey = KbName('LeftArrow');
upKey = KbName('UpArrow');
downKey = KbName('DownArrow');
if scan == 1,
    if strcmpi(os, 'mac'),
        rightKey = 33;% button box, mac: 33, windows: 52
        waitKey = 52;%trigger mac: 52 windows : 222
    elseif strcmpi(os, 'win'),
        rightKey = 52;
        waitKey = 222;
    end
elseif scan == 0,
    rightKey = KbName('RightArrow'); %54:buttonbox;
    waitKey = KbName('k');%222;%KbName('k'); %222: triggerbox;
end
[ keyIsDown, time, keyCode ] = KbCheck;

[w, rect] = Screen('OpenWindow', screenNumber, 0,[], 32, doublebuffer+1);
% [w, rect] = Screen('OpenWindow', 0, 0,[], 32, doublebuffer+1);

% alpha blending...to smooth points
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
[center(1), center(2)] = RectCenter(rect);
drawcenter=[0 0];
fps=Screen('FrameRate',w);      % frames per second
ifi=Screen('GetFlipInterval', w);
if fps==0
    fps=1/ifi;
end;

black = BlackIndex(w);
white = WhiteIndex(w);
red = [255 0 0];
HideCursor;	% Hide the mouse cursor
Priority(MaxPriority(w));

% Do initial flip...
vbl=Screen('Flip', w);

%% ----------- initialize dot positions and velocities----------------------

ppd = pi * (rect(3)-rect(1)) / atan(mon_width/v_dist/2) / 360;    % pixels per degree
pfs = dot_speed * ppd / fps;                            % dot speed (pixels/frame)
dot_size = dot_w * ppd*2;                                        % dot size (pixels)
fix_cord = [center-fix_r*ppd center+fix_r*ppd];         % center fix coord
line_width  = 3;
pt= 12;

%%------xy coord------
rmax = max_d * ppd + 150;	% maximum radius of annulus (pixels from center) %% fullfield
rmin = min_d * ppd; % minimum
r = rmax * sqrt(rand(ndots,1));	% r
r(r<rmin) = rmin;
t = 2*pi*rand(ndots,1);                     % theta polar coordinate
cs = [cos(t), sin(t)];
xy = [r r] .* cs;   % dot positions in Cartesian coordinates (pixels from center)

%----- change
ndots_one_dir = p_dir*ndots;
ndots_4dir = ndots -ndots_one_dir;
mdir = [ones(ndots_one_dir, 1); (2*rand(ndots_4dir,1)-1);]; %[real_dots; noise_ndots;]
dr_tm = pfs * mdir;                            % change in radius per frame (pixels)
dr = pfs * mdir;            % just to predefine dr
dt_tm = [ones(ndots_one_dir,1); (2*rand(ndots_4dir,1)-1);]*0.01; % [real_dots; noise_ndots];
dt = [ones(ndots_one_dir,1); (2*rand(ndots_4dir,1)-1);]*0.01; % just to predefine dt
dxdy = [dr_tm dr_tm] .* cs;                       % change in x and y per frame (pixels)
drdt_sign = [1 1; 1 -1; -1 1; -1 -1];

%% ----color vector
fname_LUT = 'isolumin/spc_mrz3t_gammaLUT.txt';
gamLUT = textread(fname_LUT);

%% ---- EDIT HERE for visual contrast and luminance -----
bg_lumi = 0.5;        % [0,1] luminance of bg
dots_contrast = .8;  %  [0,1] contrast of dots.
% ---- END EDIT --------------

% --- calculate dot-luminances randomly within contrast-range, gamma-correct them:
col_bg = img_gammaConvert(gamLUT,round(255*[bg_lumi, bg_lumi, bg_lumi])); % gamma-corrected bg color

col_high = bg_lumi*(1+dots_contrast); % cols of each set of dots.
col_low  = bg_lumi*(1-dots_contrast);

% dots should have random luminance within range: [col_low,col_high]
%       ie. the maximal contrast = dots_contrast.

colvect = col_low*255 + (col_high-col_low)*255*[rand(1, ndots);];
colvect = repmat (colvect, 3, 1);
colvect = colvect'; % [ndots, 3]
colvect = img_gammaConvert(gamLUT, round(colvect)); % gamma corrected dot-RGBs
colvect = colvect'; % [3, ndots]

%% ----dot size
if (differentsizes>0) % -Create a vector with different point sizes for each single dot, if requested:
    s=(1+rand(1, ndots)*(differentsizes-1))*dot_size;
end;

%% parameter
dx = 0;    dy = 0;
% time
% dir_change_t = [];
color_change_t = [];
resp_t = [];            dir_resp_t = [];
frame_t = [];

% number
resp_n = 0;             dir_resp_n = 0;
color_change_n = 0;     dir_change_n = 0;
color_change_one = 0;   dir_change_one = 0;
yes_resp = 0;           no_resp = 0;
scan_n = 0; scan_t = 0;
color_change_i = [];    dir_change_i = [];
color_indx_array = [];  dir_change_array = [];
color_resp_array = [];        dir_resp_array = [];
result = [];
save_on = 1;
resp = [];

nframes = attention_frame*n;
attention_switch_t = [];

%% TEXT
textSize=24;%24;
htxtshift=textSize/2; htxtmore = 20;
vtxtshift=textSize/2; vtxtmore = 10;%textSize/2;

% ------------------------------ ANIMATION LOOP---------------------------
% if (doublebuffer==1)
%     vbl=Screen('Flip', w, vbl + (waitframes-0.5)*ifi);
% end;
Screen('TextSize',w, textSize);%30);
Screen('TextStyle', w, 0);
% Screen('DrawText', w, 'X : red', center(1), center(2), [255, 0, 0, 255]);
% Screen('DrawText', w, 'X : Clockwise, Outward' , center(1), center(2) +30, red);
if (doublebuffer==1)
    vbl=Screen('Flip', w, vbl + (waitframes-0.5)*ifi);
end;

% waitsecs(1);
attention_sess = 0;
telapsed = 0;
scan_time = [];
% figure;
fix_color = [0 0 255];
%% When the program is initiated, it will show green SCREEN!!!
Screen('FillRect', w , [128 200 128] );
Screen('Flip', w);

%--- eyetracker --- %
if eyetrack ==1,
    % initialize
    pathofdll='D:\Arrington Eye Tracker\ViewPoint (2.8.4 beta 33d)\VPX_InterApp.dll' ;
    pathofheader1='D:\Arrington Eye Tracker\ViewPoint (2.8.4 beta 33d)\SDK\vptoolbox.h';
    pathofheader2='D:\Arrington Eye Tracker\ViewPoint (2.8.4 beta 33d)\SDK\vpx.h';
%     vpx_initialize(PathOfDll,PathOfHeader1, PathOfHeader2) % natalia.. was uncommented
    vpx_initialize(pathofdll,pathofheader1, pathofheader2) % natalia

    str = sprintf('%s', 'VideoMode Precision'); % 30 Hz
    str = sprintf('%s', 'VideoMode Speed'); % 60 Hz
    VPX_SendCommandString(str);

    open data file
    str = sprintf('%s', 'dataFile_Pause');
    VPX_SendCommandString(str);

    configure synchronous string insertion
    str = sprintf('%s', 'dataFile_asynchStringData No');
    VPX_SendCommandString(str);

    thedatestr=datestr(now,'yyyy-mm-dd_HH.MM');
    thestr_vpx = ['"vpx_atten_', subject, '_',sess, num2str(sess), thedatestr,'.txt"'];
    str = sprintf('dataFile_NewName %s', thestr_vpx );
    vpx_SendCommandString (str);
end
%--- eyetracker --- %
if scan
    wait_MRtrigger(waitKey);
end

tic
Screen('FillRect', w , [128 128 128] );
Screen('Flip', w);
tic
% Screen('DrawText', w, 'X', center(1), center(2), [128 128 128]);
Screen('Flip', w);
% waitsecs(10);
if eyetrack,
    str = sprintf('%s', 'dataFile_Resume');
    VPX_SendCommandString(str);
end
for i = 1:nframes
    telapsed = toc;
    Screen('FillRect', w, col_bg)
    if mod(i, attention_frame) == 1 % attention switch
        %         attention = mod(attention + 1, 3);
        attention_sess = attention_sess + 1;
        attention = sequ(attention_sess);
        attention_switch_t = [attention_switch_t; attention telapsed;];
        if eyetrack,
            str = sprintf('%s', 'dataFile_InsertString CONDITION#', num2str(attention));
            VPX_SendCommandString(str);
        end
    end

    if (i>1)
        Screen('TextSize', w, textSize);
        %         Screen('TextStyle', w, 0);
        Screen('TextFont',w, 'Times');
        if attention == 1

            if mod(i, attention_frame) < 120
                Screen('DrawText', w, 'TASK', center(1)-htxtshift-htxtmore, center(2)-vtxtshift-vtxtmore, fix_color);
            else
                if mod(i,10) <10
                    Screen('DrawText', w, 'T', center(1)-htxtshift, center(2)-vtxtshift-vtxtmore, fix_color);
                end
            end
            Screen('DrawDots', w, xymatrix, s, colvect, center,1);  % change 1 to 0 to draw square dots
        elseif attention == 0
            if mod(i, attention_frame) < 120
                Screen('DrawText', w, 'REST', center(1)-htxtshift-htxtmore, center(2)-vtxtshift-vtxtmore, fix_color);
            else
                if mod(i,10) <10
                    Screen('DrawText', w, 'X', center(1)-htxtshift, center(2)-vtxtshift-vtxtmore, fix_color);
                end
            end
            Screen('DrawDots', w, xymatrix, s, colvect, center,1);  % change 1 to 0 to draw square dots
        end
    end;

    [ keyIsDown, time, keyCode ] = KbCheck;
    if keyCode(escapeKey)
        save_on = 0;
        break;
    end;

    % coord
    cs = [cos(t), sin(t)];
    xy = [r r] .* cs;   % dot positions in Cartesian coordinates (pixels from center)

    xy = xy + dxdy;						% move dots
    t = t + dt;                         % theta change
    r = r + dr;							% update polar coordinates too

    %% --check to see which dots have gone beyond the borders of the annuli
    %- original coord
    %% ======================= FULL FIELD============================
    %     r_out = find(r > rmax | r < rmin | rand(ndots,1) < f_kill);	% dots to reposition
    r_out = find(r < rmin | rand(ndots,1) < f_kill);	% dots to reposition
    nout = length(r_out);

    if nout
        % choose new coordinates

        r(r_out) = rmax * sqrt(rand(nout,1)); %% full field
        %         r(r_out) = rmax * sqrt(rand(nout,1));
        r(r<rmin) = rmin;
        t(r_out) = 2*pi*(rand(nout,1));

        % now convert the polar coordinates to Cartesian

        cs(r_out,:) = [cos(t(r_out)), sin(t(r_out))];
        xy(r_out,:) = [r(r_out) r(r_out)] .* cs(r_out,:);

        % compute the new cartesian velocities

        dxdy(r_out,:) = [dr_tm(r_out) dr_tm(r_out)] .* cs(r_out,:);
    end;

    %% ----------dir change - to update dr, dt-----------------------------
    if (mod(i, dir_change_rate) == 1) % To vary STD in every dir_change_rate
        dir_change_variance_array = randperm(dir_change_variance);
    end
    if mod(dir_change_n,3) == 0
        dir_change_tm = randperm(4);
    end
    if mod(i, dir_change_rate) <= (dir_change_variance + 1)
        if mod(i, dir_change_rate) == (dir_change_variance_array(1) + 1)
            %             dir_indx = dir_indx(1);
            dir_indx = dir_change_tm(mod(dir_change_n, 4) + 1);
            % p_dir
            ndots_one_dir = round(p_dir*ndots);
            ndots_4dir = ndots -ndots_one_dir;
            mdir = [ones(ndots_one_dir, 1); (2*rand(ndots_4dir,1)-1);]; %[real_dots; noise_ndots;]
            dr_tm = pfs * mdir;                            % change in radius per frame (pixels)
            dr = pfs * mdir;            % just to predefine dr
            dt_tm = [ones(ndots_one_dir,1); (2*rand(ndots_4dir,1)-1);]*0.01; % [real_dots; noise_ndots];
            dt = [ones(ndots_one_dir,1); (2*rand(ndots_4dir,1)-1);]*0.01; % just to predefine dt
            dxdy = [dr_tm dr_tm] .* cs;
            dr = drdt_sign(dir_indx,1)* dr_tm;
            dt = drdt_sign(dir_indx,2)* dt_tm;
            dir_change_array = [dir_change_array; i telapsed dir_indx];

            if dir_change_n > 0
                if dir_indx == 1 & dir_indx ~= dir_change_array(end-1)% if dir == DIR1 & except when 1st dir is DIR1
                    if attention == 1
                        dir_change_i = [dir_change_i; i];
                        dir_change_t = telapsed;
                        dir_change_one = dir_change_one + 1;
                    end
                end
            end
            dir_change_n = dir_change_n + 1;
        end
    end



    %---------------Response to DIRECTION change --------------------------
    if attention == 1
        if dir_change_one > 0 % if DIRECTION is changed
            if i == dir_change_i(dir_change_one)
                check_yes = 0; check_no = 0;
            end
            if i > dir_change_i(dir_change_one) + dir_change_det_min & i < dir_change_i(dir_change_one) + dir_change_det_th % 1-59 frames
                if keyCode(rightKey)
                    check_yes = check_yes + 1;
                    if dir_resp_n == 0
                        dir_resp_t = [dir_resp_t; telapsed];
                        dir_resp_n = dir_resp_n + 1;
                    elseif dir_resp_n > 0
                        if telapsed - dir_resp_t(dir_resp_n) > 0.3
                            dir_resp_t = [dir_resp_t; telapsed];
                            dir_resp_n = dir_resp_n + 1;
                        end
                    end
                else
                    check_no = check_no + 1;
                end
            end
            if i == dir_change_i(dir_change_one) + dir_change_det_th
                if check_yes > 0
                    yes_resp = yes_resp + 1;
                    dir_resp_array = [dir_resp_array; dir_resp_t(end) 1];
                else
                    no_resp = no_resp + 1;
                    dir_resp_array = [dir_resp_array; NaN 0];
                end
                yes_dir = sum(dir_resp_array(:,2)==1)/size(dir_resp_array,1);
                result = [result; attention dir_change_t(end, :) dir_resp_array(end, :) dir_resp_array(end, 1) - dir_change_t(end, 1)]
                % Quest
                resp_dir = result(end,4);
                % apply Quest only in the first training session
                if strcmpi(sess, 'psych')
                    intensity_dir = QuestQuantile(q_dir);
                    if intensity_dir >= 0,
                        intensity_dir = -0.01;
                    end
                    q_dir = QuestUpdate(q_dir, intensity_dir, resp_dir);

                    p_dir_array = [p_dir_array; p_dir];
                    p_dir = 10^intensity_dir;
                end
            end
        end

    end
    % scan time stamp
    if  keyCode(waitKey)
        if scan_n == 0% 1st response
            scan_t = [scan_t; telapsed;];
            scan_n = scan_n + 1;
        elseif scan_n > 0 % 2nd resp - last resp
            if telapsed - scan_t(end) > 0.5
                scan_t = [scan_t; telapsed;];
                scan_n = scan_n + 1;
            end
        end
    end
    if keyCode(rightKey)
        if size(resp,1) ==0
            resp = [resp;  attention dir_indx telapsed ];
        elseif telapsed - resp(end, 3) > 0.3
            resp = [resp;  attention dir_indx telapsed ]
        end
    end
    %--------------

    % matrix & buffer
    xymatrix = transpose(xy);
    if (doublebuffer==1)
        Screen('Flip', w);
    end;

    frame_t = [frame_t; telapsed];
end;
Priority(0);
ShowCursor
Screen('CloseAll');

if save_on ==1,
    save(['../02_PsychResult/resp_vis_', subject,'_',sess,'.mat']);
    if strncmpi(sess, 'sess', 4),
        runindex = runindex + 1;
        save runindex runindex
    end
end

if eyetrack,
    VPX_SendCommandString(str);
    str = sprintf('%s', 'dataFile_Close');
    VPX_SendCommandString(str);
    vpx_unload;
end