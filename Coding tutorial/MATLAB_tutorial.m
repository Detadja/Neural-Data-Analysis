%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Bioen 460/560 / EE 460/560 / CSE 490N MATLAB Tutorial 

%In this tutorial, you will learn the basics of MATLAB and skills relevant
%to completing Project 2 (the analysis of ECoG data). 

%Written by Luke M Bun 10-16-24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Section 1. Variables and indexing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB (and Python) are basically fancy graphing calcultors. You can set
% variables, manipulate, and display them. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Variables %%%
% Variables are important! 
a = 1 %without semicolon
b = 2; %with semicolon 
you_can_use_whatever_letters_you_want = 460
or_numbers_1234 = 560
%Note: you CANNOT start variable names with numbers

% Using variables (it's just math!) 
a+b 
a*b
a/b
a^b 

%% Arrays %%%
% Variables can also be composed of multiple numbers. A variable of
% multiple values is called an array. 
a = [1,2,3,4] 
b = 1:4 %do the above automatically
c = 1:0.5:4 %with smaller spacing 
d = linspace(1,4,8) %specify the number of values; 

%Some arrays can even be multidimensional 
%Here's a 2d array
twoD_array = [1,2; 3,4]
another_array = [1,2
                 3,4]
size(twoD_array) %give me the dimensions of the array

%TODO: How would you make a 3x3 array?

%% Indexing %%%
%What if I want a particular value in my array? 
a(1) %give me the first value of a 
a(1:2) %give me the first and second values
a(end) %give me the last value 
a(end-1) %give me the second to last value
a(1:end) %give me every value between the begenning and end 
%a(0.5) %why won't this work?  

%Indexing in a multidimensional array
disp(twoD_array(1,2)) %give me the value in the first row, second column

%TODO: How can you use indexing to subtract the 1st element of a from the
%4th?

%% Logical indexing %%% 
%Sometimes you need more flexibility than normal indexing. 
%One of these methods is called "logical indexing." This method uses
%booleans, which are arrays where each value is 0 (false) or 1 (true)

%Example:
%Set-up
t = 1:10; %10 time points 
y = linspace(41,50,10); %10 samples between 41 and 50
t1 = 2; 
t2 = 8; 

%Using inequalities
Lt1 = t>t1 %this is known a boolean (or logical array)
Lt2 = t<t2 
Lt = Lt1 & Lt2 %you can use & to combine logical arrays to identify where they overlap

disp(t) %original time 
disp(Lt) %true and false 
disp(t(Lt)) %times kept based on true and false

disp(y) %original data values
disp(y(Lt)) %give me all the values of y between t1 and t2 using the logical

%% Other types of variables %%% 
string = 'You can store words or sentences'; %these are called strings
combo_string = ['you can also ', 'combine strings like this']
num2str(460) %you can turn numbers into strings 
str2double('560') %or strings into numbers

mix_n_match = {string,a,c}; %These are called cells and can store strings, arrays and other cells
disp(mix_n_match{1}) %use {} to index into a cell

struct.field = mix_n_match; %Structures let you put even more things together
struct.something = a; %each part of a structure is called a field
struct.another_thing = c; %you can give each field whatever name you want

disp(struct)
disp(struct.something)

%% Special variables %%%
%Some variables come pre-assigned by MATLAB. Highlights include: 
% - pi = 3.1415...
% - inf = infinity 
% - nan = not a number 
% - eps = the smallest number MATLAB can represent
% - i or j = complex number (sqrt(-1)) 

%%% Comments %%%
% You've probably noticed that all of my explanations start with a %. In
% MATLAB, % denote comments or parts of lines that aren't code. 
%
% Using %% at the start of a line starts a section, which is a module of
% code that can be run all together. 
%
% You can also comment out multiple lines with %{ and %} at the start of
% lines. %{ is the start and %} is the end. 
%{ 
 Multi
 line 
 comment 
%} 

%% Section 2. Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Functions are operations you can perform on variables to make new
%variables. Many come pre-installed in MATLAB but you can install more or
%write your own! We've already used some: linspace and disp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Common functions %%%
input = randn(10,1); %The rand function gives you an array of N random numbers 
                  % pulled from a Gaussian distribution with mean 0 and
                  % standard deviation 1. 
output1 = max(input); %The max function gives you the biggest number in an 
                     %array 
output2 = min(input); %min gives the smallest number 
output3 = abs(input); %the abs function takes the absolute value 
output4 = median(input); %the median gives the median 

%If you don't know how a function works, you can use the help feature! 
%Ex: type "help randn" in the Command Window
%Google is also your friend: https://www.mathworks.com/help/matlab/ref/randn.html 

%TODO: What other simple calculations and functions do you think exist? How 
%would you look them up? 
%% Plotting %%%
%Visualizing your data is the critical to understanding it and debugging.
%Plotting is just another instance of a function. 

x = -10:0.1:10; 
y = -x.^2+1; %MATLAB follows PEMDAS
%(Why did I need to do .^ and not ^? What does MATLAB say if I just did ^?) 

figure; %Set up a figure to plot in 
plot(x,y); %make the plot
title('My first parabola') %make a title
xlabel('x values') %label your axis
ylabel('y values') %label your axis

%TODO: try your own formulas and make your own plots

%% Plotting 2: Histograms %%%
%There are also other types of plots you may need. One of them is histograms. 
%Histograms represent the number of times a vaue appears in an array. Let's use 
%a histogram to see if randn produces a normal distribution (bell curve). 

% Generate random numbers from a normal distribution
x = randn(1000, 1); %What happens to the distribution as you use more random numbers?

% Create the first subplot with default bins
figure; %make a new figure to plot stuff on
subplot(1, 2, 1); % create a 1x2 grid of subplots and select the first
histogram(x, 100); % you can specify the number of bins
ylabel('Count');
xlabel('Value');
title('Histogram of randn');

% Create the second subplot with custom bins
custom_bins = -4:0.5:4; % use 0.5 spacing
subplot(1, 2, 2); % select the second subplot
histogram(x, custom_bins); % you can also specify the bin edges
ylabel('Count');
xlabel('Value');
title('Histogram of randn');

% Lastly, you can also get the counts of each bin
counts = histcounts(x, custom_bins);
disp(counts);

%% Other useful functions %%%
%For Project 2, you may also want to use the findpeaks command. This does
%not come preinstalled on all versions of MATLAB, so you may need to
%install the "signal toolbox." Check if you have it by typing "which
%findpeaks" in the Command indow. 
%
% If you don't, go to the top bar, go to home and click on Add-Ons. Look up
% the "signal toolbox" and download it. You should have access to the
% function afterwords.

%Findpeaks 
x = linspace(0,6*pi); %the pi variable is automatically set to the value of pi
y = sin(x); 
[peaks, peak_indicies] = findpeaks(y); %findpeaks can have multiple outputs
disp(peaks)
disp(peak_indicies)

figure; 
plot(x,y)
hold on %lets you plot multiple things on top of each other
plot(x(peak_indicies),peaks,'x','linewidth',2)
hold off
xlabel('x')
ylabel('y')
title('findpeaks in action')

%findpeaks also has many other useful settings (like setting thresholds) so
%be sure to look up how these functions work!!!

%% Custom functions %%%
%You can also make your own functions. In MATLAB, functions need to go at
%the end of the script, so move forward to Section 6. Here's an example with
%a custom ReLU (rectified linear unit) function
x1 = -10:10; 
[x2,~] = relu(x1'); %If you don't want an output, you can replace the variable name with ~
disp(x1); 
disp(x2); 

figure; 
plot(x1,x2)
xlabel('x1')
ylabel('x2')
title('relu')

%TODO: make your own function and use it

%% Using the filter500to5k function %%%
%Custom Functions can also be saved as their own file and used. The
%provided filter500to5k function is one such example. 
%Note: to use it, you should rename it to not have the "filter500to5k"
%part. 

%Practice data 
t = linspace(0,6*pi,20000); 
test_data = sin(t) + randn(1,20000)/4; 
fs = 20000; %sampling frequency in Hz

%Change your directory to where you put the function 
cd('C:\Users\lukebun\Documents\Graduate\Bioen 460') %CHANGE THIS TO WHERE YOU PUT THE FUNCTION 
[SOS,G] = filter500to5k(fs); %run the function to get its outputs
filt_test_data = filtfilt(SOS,G,test_data); %use the outputs to filter your data

%Plot it!
figure; 
plot(test_data)
hold on %lets you plot multiple things on top of each other
plot(filt_test_data)

hold off 

%% Section 3. Loops and if statements %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%You will often need to perform the same series of operations multiple
%times. Instead of writing the same code multiple times, use a loop! 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% For loops %%%
%For loops let you iterate over an array to perform a repetative action
iters = 1:10; 
store = zeros(1,10); %an array of all 0s where we will store data
for i = iters %i is a value that is changed on each loop to be the next value in iters
    store(i) = i^2; %square each i and store it away
end %end indicates where the loop ends
disp(store)

%TODO: write a for loop that multiplies the ith element of an array by the 
%ith+1 element. 

%% If statements %%%
%You likely won't need these for project 2, but they're nice to know about
%If statements let you do something when certain conditions are met. 
vals = [-1,0,1]; 
if sum(vals) < 10 
    disp('Sum of vals is less than 10')
end %if statements also need ends
if sum(vals) == 10 
    disp('Sum of values is equal to 10')
end 

%% Nested for loops and if statements %%% 
%You can have nested for loops and if statements to pull off complex
%operations. In this example, we'll loop through every value in a 2D array 
%and store values that are odd 

x = [1,2;3,4]; 
odd_vals = zeros(2,2); 
for i = 1:2 %for each row

    for j = 1:2 %for each column
        
        val = x(i,j); 
        if mod(val,2) %checks if a val is odd using the remainder after dividing by 2
            odd_vals(i,j) = val; 
        end %end if statements

    end %end inner loop

end %end outer loop
disp(odd_vals)

%Pro-tip: when nesting loops and if statements, you should label where each
%one ends for your own sanity

%% Section 4. Starting Project 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Below is some code to get you started on Project 2. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loading data %%% 
path = 'C:\Users\lukebun\Documents\Graduate\Bioen 460'; %CHANGE THIS TO YOUR DATA FOLDER!!!
cd(path); %change the directory; 
filename = 'ECoGData_AvailabletoStudents'; %name of the data
ecog = load(filename); %load the data

disp(ecog.data) %look at the data field 
fs = ecog.data.SamplingFreq; %sampling frequency 
sig = ecog.data.Signal; %cell of the signal 
stim_t = ecog.data.StimTimes; %cell of stim times

disp(['Size of signal cell: ', num2str(size(sig))])
disp(['Size of stim_t cell: ', num2str(size(stim_t))])

% General tips for Project 2: 
% - Google and the help feature are your friends! Even the most experienced
% programmers need to look up how functions and operations work. 
% - Write out the things you want your code step by step (pseudo-code). You 
% don't need to use the right commands or worry about indexing, just think
% through how you want your code to work. 
% - Come to office hours. When you run into a problem coding, it can be
% difficult to debug over email. The best time to get help is during office
% hours. 

%% Section 5. Another example %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This example is for you to look over and use as inspiration for Project
% 2. Some parts will be similar, others will be up for you to do on your
% own. 

% I collected these data from a rhesus macaque being presented
% with visual grating stimuli at different orientations (Example:
% https://www.djmannion.net/psych_programming/_images/draw_gratings_6.png).
% At the same time, I was using a single unit tungsten electrode to record
% the timing of a neuron's action potentials (spikes). 

% This code will load those data and analyze them to determine the
% preferred orientation of the neuron I recorded from.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data %%%
cd('C:\Users\lukebun\Documents\Graduate\Bioen 460') %CHANGE THIS TO WHERE YOU PUT THE DATA!!
load('orientation_data.mat')
disp(ex_data)
spikes = ex_data.sig; 
stim_on = ex_data.stim_on; 
stim_off = ex_data.stim_off; 
orient_trial = ex_data.orient; 

%Each row is a new trial where a new stimulus is presented
n_trials = size(stim_on,1); %use the size of the array to get the size 

%% Make a raster plot %%% 
figure; 
hold on 
for iT = 1:size(stim_on,1) %for each trial
    spike_t = spikes{iT,1} - stim_on(iT); %get spike times in a trial, relative to stimulus onset 

    %Plot each spike as a | at a given time along the same row
    plot(spike_t, iT*ones(length(spike_t),1),'k|') %the k specifies the color as black, | is the symbol
end 
xlim([-0.5,1.5]) %Limit the x-axis
xlabel('Time')
ylabel('Trial #')
title('Raster plot')
%What does t=0 represent? 
%What happens to the spike rate at t=0? 

%% Get spike rate for each trial %%%
spike_rates = zeros(n_trials,1); 
for iT = 1:size(stim_on,1) %for each trial
    spike_t = spikes{iT,1}; 

    %only keep the spikes from when the stimulus was on
    spike_t = spike_t(spike_t > stim_on(iT) & spike_t < stim_off(iT)); %logical indexing

    n_spikes = length(spike_t); %number of spikes I kept
    spike_rates(iT) = n_spikes / (stim_off(iT) - stim_on(iT)); %n spikes/time
end 
figure; 
plot(spike_rates)
xlabel('Trial #')
ylabel('Spike rate (sp/s)')
title('Spike rate on each trial')

%% Get mean spike rate for each orientation %%%
orients = unique(orient_trial); %find all the unique orientations presented
mean_spike_rates = zeros(length(orients),1); 
for iO = 1:length(orients) %for each orientation
    %get the spike rates on every trial with that particular orientation
    ori_spike_rates = spike_rates(orient_trial == orients(iO)); 
    mean_spike_rates(iO) = mean(ori_spike_rates); %take the mean
end
figure; 
plot(orients,mean_spike_rates)
xlabel('Stimulus orientation')
ylabel('Spike rate (sp/s)')
title('Orientation tuning curve') %A plot of neuronal activity vs stimulus property is also known as a tuning curve
%Challenge: how would you represent the standard deviation for each
%stimulus condition? (hint: use errorbar) 

%% What is the preferred orientation? %%%
%You can use the second output of the max function
[~,pref_ori_index] = max(mean_spike_rates); %second output is the index of the max value 
disp(['The preferred orienation is ', num2str(orients(pref_ori_index))]); 

%But you'll notice that there are actually two peaks: one big and one
%small. What if we wanted both of them?  
%You can try findpeaks!
[~,pref_ori_indices] = findpeaks(mean_spike_rates); 
disp(['The first peak is at ', num2str(orients(pref_ori_indices(1)))]); 
disp(['The second peak is at ', num2str(orients(pref_ori_indices(2)))]); 

%% Section 6. Custom functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions need to go at the bottom of MATLAB scripts, so this section is 
% out of order from Section 2. You can also make functions their own scripts
% (like the provided filter500to5k function). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example function %%%
function [x2,x1] = relu(x1) %one function can have multiple outputs
    %This function produces the relu operation, which sets all negative
    %numbers to 0 and returns all positive numbers unchanged. 

    x2 = x1; %make a copy of x1
    x2(x1<0) = 0; %using logical indexing to change values  
end 

