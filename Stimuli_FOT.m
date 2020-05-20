 clear all;
close all;  
 
if usejava('System.time')
    disfp('error loading java package "System.time"');
    return; 
end
if usejava('java.util .LinkedList')   
    disp('error loading');  
    return; 
end
    %x for subnum and counterbalance; create book number 0 doesn't exist... look into the correct way to format three input boxes
   
    prompt = {'Subnum' 'counterbalance'};
    name = 'q';
    numlines = 1;
    defaultanswer = {'0','1'};
    answer = inputdlg(prompt, name, numlines, defaultanswer);
    if(size(answer) ~= 2) 
        clear;
        clc;
         disp('Exiting.'); % program exits here because of booknum...
        return;
    end
    
     Priority(2);
    [subject,counterbalance] = deal(answer{:});  
    counterbalance = str2num(counterbalance);
    cbrange = (1:2);
    if counterbalance < 1 || counterbalance > 2 
        msg = 'An incorrect counterbalance was entered. Please enter a 1 or 2.';
        error(msg)
    end
PsychDefaultSetup(1);
condition = 1;
% which order?
orderArray = {'1a','1b','1c','1d','1e','1f','2a','2b','2c','2d','2e','2f','3a','3b','3c','3d','3e','3f'}; 
[selectionIndex3, leftBlank] = listdlg('PromptString', 'Select an order to run:', 'SelectionMode', 'single', 'ListString', orderArray);
order= orderArray{selectionIndex3};
%setup screen info and various associated types of information
whichScreen = 2; %allow to choose the display if there?s more than one
white = WhiteIndex(whichScreen); %pixel value for white
black = BlackIndex(whichScreen); %pixel value for black
gray = GrayIndex(whichScreen); %pixel value for middle gray
[wPtr, myRect] = Screen( 'OpenWindow', whichScreen, gray);
Screen('BlendFunction',wPtr ,GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

%set up the fixation cross, the cross is actually drawn within the task loop
[xCenter, yCenter] = RectCenter(myRect);
fixCrossDimPix = 60;
% set the fixation and center coordinates 
xCoords = [-fixCrossDimPix fixCrossDimPix 0 0];
yCoords = [0 0 -fixCrossDimPix fixCrossDimPix];
allCoords = [xCoords; yCoords];

% Set the line width for our fixation cross
lineWidthPix = 6;

%the below information is set to help with smoothing the moving image 
%Sets the units of your root object (screen) to pixels
set(0,'units','pixels');
% %Obtains this pixel information
Pix_SS = get(0,'screensize');
% % use the above and below info for smooth moving image
moveDist = Pix_SS(4)/2.2; %the value 4 represents the height in pixels, the divisor may be manipulated based on the height of the screen 
%the above statement will determine how far down the screen the image moves
refresh = Screen('GetFlipInterval', wPtr); %screen refresh rate
refresh2 = refresh *2; %cuts the speed of the image in half
vbl = Screen('Flip', wPtr);



%chooses image resource folder based on counterbalance
if strcmp(order, '1a') 
    picLoc = 'C:\Users\vpixx\Desktop\Floating Faces Task\faces folder\resized/1a';
elseif strcmp(order, '1b')
    picLoc = 'C:\Users\vpixx\Desktop\Floating Faces Task\faces folder\resized/1b';
elseif strcmp(order, '1c')
    picLoc = 'C:\Users\vpixx\Desktop\Floating Faces Task\faces folder\resized/1c';
elseif strcmp(order, '1d')
    picLoc = 'C:\Users\vpixx\Desktop\Floating Faces Task\faces folder\resized/1d';
elseif strcmp(order, '1e')
    picLoc = 'C:\Users\vpixx\Desktop\Floating Faces Task\faces folder\resized/1e';
elseif strcmp(order, '1f')
    picLoc = 'C:\Users\vpixx\Desktop\Floating Faces Task\faces folder\resized/1f';
elseif strcmp(order, '2a') 
    picLoc = 'C:\Users\vpixx\Desktop\Floating Faces Task\faces folder\resized/2a';
elseif strcmp(order, '2b')
    picLoc = 'C:\Users\vpixx\Desktop\Floating Faces Task\faces folder\resized/2b';
elseif strcmp(order, '2c')
    picLoc = 'C:\Users\vpixx\Desktop\Floating Faces Task\faces folder\resized/2c';
elseif strcmp(order, '2d')
    picLoc = 'C:\Users\vpixx\Desktop\Floating Faces Task\faces folder\resized/2d';
elseif strcmp(order, '2e')
    picLoc = 'C:\Users\vpixx\Desktop\Floating Faces Task\faces folder\resized/2e';
elseif strcmp(order, '2f')
    picLoc = 'C:\Users\vpixx\Desktop\Floating Faces Task\faces folder\resized/2f';
elseif strcmp(order, '3a') 
    picLoc = 'C:\Users\vpixx\Desktop\Floating Faces Task\faces folder\resized/3a';
elseif strcmp(order, '3b')
    picLoc = 'C:\Users\vpixx\Desktop\Floating Faces Task\faces folder\resized/3b';
elseif strcmp(order, '3c')
    picLoc = 'C:\Users\vpixx\Desktop\Floating Faces Task\faces folder\resized/3c';
elseif strcmp(order, '3d')
    picLoc = 'C:\Users\vpixx\Desktop\Floating Faces Task\faces folder\resized/3d';
elseif strcmp(order, '3e')
    picLoc = 'C:\Users\vpixx\Desktop\Floating Faces Task\faces folder\resized/3e';
elseif strcmp(order, '3f')
    picLoc = 'C:\Users\vpixx\Desktop\Floating Faces Task\faces folder\resized/3f';
else
    msg = 'Please choose a valid order number';
    error(msg)
end
%below are the imread and alpha functions to get rid of black background of
%images, load the image data into a matrix and draw a texture of that image

% [imMatrix_it1, ~, imdata_alpha] = imread(fullfile( picLoc , 'it1.png' ));
% imMatrix_it1(:,:,4)=imdata_alpha; % added alpha layer to 2 because images are greyscale
% imTexture_it1 = Screen('MakeTexture',wPtr, imMatrix_it1);
% I took these it images out as requested by Dr. Scott
% [imMatrix_it2, ~, imdata_alpha] = imread(fullfile( picLoc , 'it2.png' ));
% imMatrix_it2(:,:,4)=imdata_alpha; % added alpha layer to 2 because images are greyscale
% imTexture_it2 = Screen('MakeTexture',wPtr, imMatrix_it2);

[imMatrix_iu1, ~, imdata_alpha] = imread(fullfile( picLoc , 'iu1.png' ));
imMatrix_iu1(:,:,4)=imdata_alpha; % added alpha layer to 2 because images are greyscale
imTexture_iu1 = Screen('MakeTexture',wPtr, imMatrix_iu1);

[imMatrix_iu2, ~, imdata_alpha] = imread(fullfile( picLoc , 'iu2.png' ));
imMatrix_iu2(:,:,4)=imdata_alpha; % added alpha layer to 2 because images are greyscale
imTexture_iu2 = Screen('MakeTexture',wPtr, imMatrix_iu2);

% [imMatrix_ct1, ~, imdata_alpha] = imread(fullfile( picLoc , 'ct1.png' ));
% imMatrix_ct1(:,:,4)=imdata_alpha; % added alpha layer to 2 because images are greyscale
% imTexture_ct1 = Screen('MakeTexture',wPtr, imMatrix_ct1);
% I took these ct images out as requested by Dr. Scott
% [imMatrix_ct2, ~, imdata_alpha] = imread(fullfile( picLoc , 'ct2.png' ));
% imMatrix_ct2(:,:,4)=imdata_alpha; % added alpha layer to 2 because images are greyscale
% imTexture_ct2 = Screen('MakeTexture',wPtr, imMatrix_ct2);

[imMatrix_cu1, ~, imdata_alpha] = imread(fullfile( picLoc , 'cu1.png' ));
imMatrix_cu1(:,:,4)=imdata_alpha; % added alpha layer to 2 because images are greyscale
imTexture_cu1 = Screen('MakeTexture',wPtr, imMatrix_cu1);

[imMatrix_cu2, ~, imdata_alpha] = imread(fullfile( picLoc , 'cu2.png' ));
imMatrix_cu2(:,:,4)=imdata_alpha; % added alpha layer to 2 because images are greyscale
imTexture_cu2 = Screen('MakeTexture',wPtr, imMatrix_cu2);

[imMatrix_u1, ~, imdata_alpha] = imread(fullfile( picLoc , 'u1.png' ));
imMatrix_u1(:,:,4)=imdata_alpha; % added alpha layer to 2 because images are greyscale
imTexture_u1 = Screen('MakeTexture',wPtr, imMatrix_u1);

[imMatrix_u2, ~, imdata_alpha] = imread(fullfile( picLoc , 'u2.png' ));
imMatrix_u2(:,:,4)=imdata_alpha; % added alpha layer to 2 because images are greyscale
imTexture_u2 = Screen('MakeTexture',wPtr, imMatrix_u2);

[imMatrix_f1, ~, imdata_alpha] = imread(fullfile( picLoc , 'f1.png' ));
imMatrix_f1(:,:,4)=imdata_alpha; % added alpha layer to 2 because images are greyscale
imTexture_f1 = Screen('MakeTexture',wPtr, imMatrix_f1);

[imMatrix_f2, ~, imdata_alpha] = imread(fullfile( picLoc , 'f2.png' ));
imMatrix_f2(:,:,4)=imdata_alpha; % added alpha layer to 2 because images are greyscale
imTexture_f2 = Screen('MakeTexture',wPtr, imMatrix_f2);


%set individual arrays for each species of faces
%imTexture_it = [imTexture_it1 imTexture_it2];
imTexture_iu = [imTexture_iu1 imTexture_iu2];
%imTexture_ct = [imTexture_ct1 imTexture_ct2];
imTexture_cu = [imTexture_cu1 imTexture_cu2];
imTexture_u = [imTexture_u1 imTexture_u2];
imTexture_f = [imTexture_f1 imTexture_f2]; 


    % Connect to NetStation
    DAC_IP = '10.10.10.42';
    NetStation('Connect', DAC_IP, 55513);
    NetStation('Synchronize');
    NetStation('StartRecording');
    
'Press any button to start'    
%set task info and fix lines
numRuns = 4;
if condition == 1
    numTrial = 4;
    nSpecies = 4;
    speciesName = 'fix';
          label = 'fx07';
   KbWait();       
   WaitSecs(3); %sets a 3 second delay before fixation appears on the screen
  
    %draw fix
    Screen('DrawLines', wPtr, allCoords,lineWidthPix, black, [xCenter yCenter], 2);   
    % Flip to the screen
    Screen('Flip', wPtr);
    speciesName %label
    tic

    NetStation('Event',label, GetSecs, GetSecs+cputime, 'trl#',condition,'race',speciesName); % signals the beginning of a fixation cross

    WaitSecs(5);     % Wait for 5 seconds to move to trial
    %the below values will be the initial values for indexing face
    %textures from the major texture array, these values will be changed in the
    %following loop to reflect the new locations of the next species' faces
    w = 1;
    x = 2;
    toc
    for trial=1:numTrial %start labeling based on the correct counterbalance
        %this imTextures array is a 1x12 array containing the above textures
        if counterbalance == 1
            imTextures = [imTexture_iu imTexture_cu imTexture_u imTexture_f]; %imTexture_it and imTexture_ct taken out of here
            for stim = 1:nSpecies 
%                 if trial == 1
%                     speciesName = 'it';
                if trial ==1
                    speciesName = 'iu';
%                 elseif trial ==3
%                     speciesName = 'ct';
                elseif trial ==2
                    speciesName = 'cu';
                elseif trial ==3
                    speciesName = 'un';
                elseif trial ==4
                    speciesName = 'fa';
                end
                    switch speciesName %Change the event label based on species
%                         case 'it'
%                             label = 'it01';
                        case 'iu'
                            label = 'iu02';
%                         case 'ct'
%                             label = 'ct03';
                        case 'cu'
                            label = 'cu04';
                        case 'un'
                            label = 'un05';
                        case 'fa'
                            label = 'fa06';
                    end
            end
        elseif counterbalance == 2
            imTextures = [imTexture_cu imTexture_iu imTexture_u imTexture_f];  %imTexture_it and imTexture_ct taken out of here
            for stim = 1:nSpecies 
                if trial == 1
                    speciesName = 'cu';
%                 elseif trial ==2
%                     speciesName = 'it';
                elseif trial ==2
                    speciesName = 'iu';
%                 elseif trial ==4
%                     speciesName = 'ct';
                elseif trial ==3
                    speciesName = 'un';
                elseif trial ==4
                    speciesName = 'fa';
                end
                    switch speciesName %Change the event label based on species
%                         case 'it'
%                             label = 'it01';
                        case 'iu'
                            label = 'iu02';
%                         case 'ct'
%                             label = 'ct03';
                        case 'cu'
                            label = 'cu04';
                        case 'un'
                            label = 'un05';
                        case 'fa'
                            label = 'fa06';
                    end
            end
        else
            msg = 'An incorrect counterbalance was entered. Please enter a 1 or 2.';
            error(msg)
        end
  
        tic %begins the timing sequence for species (send this out to a document eventually)
        numFaces = [w x]; %set 2 faces for this species 1 and 2 in the main array
        for N = numFaces(1):numFaces(2) % for 1st face through the 2nd face, present the following code, iterating through imTextures array to grab new faces
        speciesName %read out species name for assurance
        NetStation('Event',label, GetSecs, GetSecs+cputime, 'trl#',trial,'race',speciesName);
        tic %begins the timing sequence for species (send this out to a document eventually)    
        numTimesFace = 1;
            for h = 1:numTimesFace
                r1=[0 0 450 500]; %draw first rect (change based off of size of image)
                r2=OffsetRect(r1, xCenter-225 ,-500);%draw 2nd rect (change the xCenter subtractor by half the size of the images)
                Screen('DrawTexture', wPtr, imTextures(N), r1, r2); %Fill the buffer with the first texture
                vbl = Screen('Flip', wPtr,  vbl + (1.5) * refresh2); %update the display with the buffer content
                WaitSecs(.007);%relative to the refresh rate
                %    %begin using loop concept to draw rectangles and draw textures to
                %    rectangles
                r_prev = r2;
                for i = 3:moveDist
                    rnew = OffsetRect(r_prev, 0,4); %new rectangle value to be used to offset rect
                    Screen('DrawTexture', wPtr, imTextures(N), [], rnew);%draw text to the new offset rect(using selected imTexture)
                    vbl = Screen('Flip', wPtr, vbl + (1.5) * refresh2);
                    r_prev = rnew; % set previous r value to the most recent rnew value, will be used for next offset rect
                    WaitSecs(.007);%wait secs is relative to the refresh rate.
                end
            end
            trialtime = toc %ends timing sequence for species
        end
        
        %change values of w, x to index face textures from the correct locations in main imTextures array
        w = w + 2;
        x = x + 2;


    %draw fixation cross
    speciesName = 'fix';     
    label = 'fx07';  %sets label to be used for fixation cross
    Screen('DrawLines', wPtr, allCoords,lineWidthPix, black, [xCenter yCenter], 2);   
    % Flip to the screen
    Screen('Flip', wPtr);
    speciesName
    tic
    % Wait for a 10 seconds to move to trial
    NetStation('Event',label, GetSecs, GetSecs+cputime, 'trl#',condition,'race',speciesName); % signals the beginning of a fixation cross
    WaitSecs(5); %cut fixaton to 5 seconds long 
    toc
    end

    %disconnect from netstation
    NetStation('Synchronize');
    NetStation('StopRecording');
    NetStation('Disconnect', '10.10.10.42');
end

%below will display the 'Thank you for participating' text
ifi = Screen('GetFlipInterval', wPtr);
[screenXpixels, screenYpixels] = Screen('WindowSize', wPtr);
Screen('TextSize', wPtr, 45);
DrawFormattedText(wPtr, 'Thank you for participating', screenXpixels * 0.4, screenYpixels * 0.5, [1 0 0]);
Screen('Flip', wPtr);
speciesName
tic
WaitSecs(5);
toc

sca;
close all;
clear all;    