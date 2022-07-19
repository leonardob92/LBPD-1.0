#%% Temporal pattern recognition in the human brain: a dual simultaneous processing - LEONARDO BONETTI

#%%


#BEFORE PROCEEDING, PLEASE NOTE:

#As follows, every analysis that we made has been reported to ensure full disclosure.
#Please, note that in this Leading script, I use built-in functions,
#functions from well-known toolboxes (e.g. OSL, SPM, FieldTrip) and
#in-house-built functions, which are reported together with this script in the collection named LBPD.
#Data can be provided according to the Danish regulation, upon reasonable request.
#If interested, please contact Leonardo Bonetti.
#More information is provided in the ReadMe.mat file that I strongly advise to read.
#Finally, in this script the memorized (M) and novel (N) patterns have been
#sometimes called old (corresponding to M) and new (corresponding to N).


#Leonardo Bonetti
#leonardo.bonetti@psych.ox.ac.uk
#leonardo.bonetti@clin.au.dk






#%% -*- coding: utf-8 -*-
"""
Curve fitting for "Temporal pattern recognition in the human brain: a dual simultaneous processing"

Leonardo Bonetti

"""


#%% IMPORTING SOME PACKAGES


import scipy.io #to read matlab data
import matplotlib.pyplot as plt
# from math import e,sin
from scipy.optimize import curve_fit
import numpy
from numpy import savetxt
# import math



#%%

# BLOCKS 1 = 0.1-1 Hz OLD; 2 = 2-8 Hz OLD; 3 = 2-8 Hz NEW


#%% GAUSSIAN SLOW NEGATIVITY

#in def gaus you have one line for skewed (log) gaussian and one line for standard gaussian

#BLOCK 1: 0.1-1 Hz, GAUSSIAN


#create csv file
filename = '/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/Tables/Table_4/Old_01_1Hz_Coeffs_R.csv'
csvfl = open(filename,'w')
csvfl.write('Parcel;R_squared;a1;b1;c1;a2;b2;c2;a3;b3;c3 \n')

#reading matlab data
mat = scipy.io.loadmat('/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/K_means/Slowneg_1_Average_l_0_Old_l_1/Clusters.mat')
#reading time
mat2 = scipy.io.loadmat('/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/time.mat')
timeti = mat2["time"]
timeseries = mat["KJK"] #extracting timeseries
no_parcels = [0,6,9,16] #parcels that are not really meaningful here..
#saving predicted data for plotting in Matlab (opening the file here)
dir2 = '/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/Tables/' #'.csv'

#actual computation
for timeseries_index in range(len(timeseries)): #over time-series (previously manually inserted index of time-series to be used from timeseries variable)
    print(timeseries_index)    
    if timeseries_index not in no_parcels:
        alltime = [1,3,4,5,8,12,13,18,11,14,7,15,17] #timeseries used with all time-points (to fit the equation)
        gaussnorm = [1,3,4,5,8,12,13,18,11,14,2,10] #timeseries with normal Gaussian function
        p0_1 = [1,3,4,5,8,12,13,18,2] #timeseries with some initial values used by the fitting function
        p0_2 = [15,17] #timeseries with some other initial values used by the fitting function
        #defining time-points to be used
        if timeseries_index in alltime:
            time_point_start = 1
            time_point_end = 390
        elif timeseries_index == 2:
            time_point_start = 1
            time_point_end = 130
        elif timeseries_index == 10:
            time_point_start = 200
            time_point_end = 390
        #defining whether using normal Gaussian function or Gaussian pulse (summation of several Gaussian functions, here 3)
        if timeseries_index in gaussnorm:
            guassian_l = 1 # 1 for single Gaussian; 0 for Gaussian pulse (repeated Gaussians) 
        else:
            guassian_l = 0 # 1 for single Gaussian; 0 for Gaussian pulse (repeated Gaussians) 
        #defining starting points for fitting, depending on the time-series used
        if timeseries_index in p0_1:
            P0 =[1,50,0.8]
        elif timeseries_index in p0_2:
            P0 = [1,2,0.8,1,250,0.8,1,300,0.8]
        elif timeseries_index == 11:
            P0 = [1,200,0.8]
        elif timeseries_index == 14 or timeseries_index == 10:
            P0 = [1,300,0.8]
        elif timeseries_index == 7:
            P0 = [1,2,0.8,1,30,0.8,1,300,0.8]
    
        #preparing time and data for fitting
        jes = list(range(time_point_start,time_point_end))
        data = timeseries[timeseries_index][jes]
        # time = list(range(1,jes[-1])) #array with time-points
        timet2 = timeti[0][jes] #time for time-series used to fit the function (so only selected time-points of that time-series)
        timet = timeti[0][0:390] #time for plotting "full" time-series
        time = jes
        jes2 = list(range(1,391)) #getting again all time-points (slightly elaborated way which could probably be implemented in less passages, but it is ok..)
        
        if guassian_l == 1:
            #NORMAL GAUSSIAN
            def gaus(x,a,x0,sigma):
                # return a*numpy.exp(-(numpy.log(x)-x0)**2/(2*sigma**2)) #skewed gaussian (obtained by log(x) instead of simply x) THIS WORKS ONLY FOR APPARENTLY NON-NEGATIVE VALUES TO BE FITTED (SO FOR INSTANCE WITH FIRST 200 TIME-POINTS IN CASE OF TIMESERIES 2)
                return a*numpy.exp(-(x-x0)**2/(2*sigma**2)) #original standard gaussian
            popt,pcov = curve_fit(gaus,time,data,P0) #timeseries 1,3,4,5,8,12,13,18; 7,15,17(several Gaussian pulses may help); #timeseries 2 (ONLY WITH TIME-POINTS 1-130)
            # popt,pcov = curve_fit(gaus,time,data,p0=[1,200,0.8]) #timeseries 11
            # popt,pcov = curve_fit(gaus,time,data,p0=[1,300,0.8]) #timeseries 14; 10(ONLY WITH TIME-POINTS 200-390)
            a,b,c = popt
            #getting values predicted by the model        
            y_line = gaus(jes2,a,b,c) #getting predicted values
            data = timeseries[timeseries_index][jes2] #getting data of full time-series (for later plotting purposes)
        else:
            #GAUSSIAN PULSE
            def gauspulse(x,a,x0,sigma,a2,x02,sigma2,a3,x03,sigma3):
                # return a*numpy.exp(-(numpy.log(x)-x0)**2/(2*sigma**2)) #skewed gaussian (obtained by log(x) instead of simply x) THIS WORKS ONLY FOR APPARENTLY NON-NEGATIVE VALUES TO BE FITTED (SO FOR INSTANCE WITH FIRST 200 TIME-POINTS IN CASE OF TIMESERIES 2)
                return a*numpy.exp(-(x-x0)**2/(2*sigma**2)) + a2*numpy.exp(-(x-x02)**2/(2*sigma2**2)) + a3*numpy.exp(-(x-x03)**2/(2*sigma3**2))
            popt,pcov = curve_fit(gauspulse,time,data,P0) #timeserie 7
            # popt,pcov = curve_fit(gauspulse,time,data,p0=[1,2,0.8,1,250,0.8,1,300,0.8]) #timeserie 15,17
            #GAUSSIAN PULSE
            a,x0,sigma,a2,x02,sigma2,a3,x03,sigma3 = popt
            #getting values predicted by the model
            y_line = gauspulse(jes2,a,x0,sigma,a2,x02,sigma2,a3,x03,sigma3) #getting predicted values
            data = timeseries[timeseries_index][jes2] #getting data of full time-series (for later plotting purposes)
        
        #getting R-squared value for model evaluation
        #residual sum of squares
        residuals = data - y_line 
        ss_res = numpy.sum(residuals**2)
        #total sum of squares
        ss_tot = numpy.sum((data-numpy.mean(data))**2)
        #R-squared value
        r_squared = 1 - (ss_res / ss_tot)
        
        #saving coefficients and R-squared in csv file (this could be done better, but it works fine and it is just slightly less elegant than other solutions..)
        if guassian_l == 1:
            lrow = '{};{};{};{};{} \n'
            lrow = lrow.format(timeseries_index + 1,r_squared,a,b,c)
            csvfl.write(lrow)
        else:
            lrow = '{};{};{};{};{};{};{};{};{};{};{} \n'
            lrow = lrow.format(timeseries_index + 1,r_squared,a,x0,sigma,a2,x02,sigma2,a3,x03,sigma3)
            csvfl.write(lrow)
            
        #plotting original data, subset of the data used for fitting function (if that is not the same as the full data), fitted function
        x_line = numpy.arange(min(timet),max(timet),1)
        # plt.plot(timet,timeseries[timeseries_index][0:390],color='orange') #plotting full timeseries
        plt.plot(timet,y_line,'--',color='red') #plotting fitted function (only within selected time-points)
        plt.plot(timet,data,color='blue') #plotting timeseries (only selected time-points)
        plt.ylim(-2.9,3.71) #y-axis limits
        plt.grid(b=True, which='both',linewidth=0.5) #grid
        plt.minorticks_on() #necessary to get minor grid
        #save figure
        pathfig = '/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/Figures/Figure_3/Timeseries/SlowNeg_Old_parcel_' + str(timeseries_index + 1)
        plt.savefig(pathfig  + '.eps',format='eps')
        plt.savefig(pathfig  + '.png',format='png')
        plt.show()
        #saving timeseries for plotting them in Matlab (to get uniform graphics)
        savetxt(dir2 + 'slowfreq_predicted_data_' + str(timeseries_index + 1) + '.txt', y_line, delimiter = ',') #predicted data (y_line)
        savetxt(dir2 + 'slowfreq_real_data_' + str(timeseries_index + 1) + '.txt', data, delimiter = ',') #real data (data)
    else:
        lrow = '{} \n'
        lrow = lrow.format(timeseries_index + 1)
        csvfl.write(lrow)
        
#closing csv file
csvfl.close()

#%%






#%% SINGLE TONE - N100 - FIVE TONES WITH MORE PROPER SUMMATION OF DIFFERENT SWEKED GAUSSIAN MULTIPLIED BY COSINE (SINUSOIDAL FUNCTION)

#%% FOR N100

#BLOCK 2: 2-8 Hz, WAVELET SUMMATION - OLD
#BLOCK 3: 2-8 Hz, WAVELET SUMMATION - NEW


#create csv file
filename = '/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/Tables/Table_6/OldNew_2_8Hz_Coeffs_R.csv'
csvfl = open(filename,'w')
csvfl.write('Parcel;Condition;R_squared;a1;b1;c1;omega1;phi1;a2;b2;c2;omega2;phi2;a3;b3;c3;omega3;phi3;a4;b4;c4;omega4;phi4;a5;b5;c5;omega5;phi5 \n')

#loading time
mat2 = scipy.io.loadmat('/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/time.mat')
timeti = mat2["time"]
timet2 = timeti[0][jes] #time for time-series used to fit the function (so only selected time-points of that time-series)
timet = timeti[0][0:315] #time for plotting "full" time-series
jes2 = list(range(1,316)) #getting again all time-points (slightly elaborated way which could probably be implemented in less passages, but it is ok..)
#defining time-points to be used for fitting purposes
time_point_start = 1
time_point_end = 315

#loading useful set of starting point for this particular equation (variable popo)
popo = numpy.load('/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/ToPython_2_8Hz/popo.npy')
#saving predicted data for plotting in Matlab (opening the file here)
dir2 = '/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/Tables/' #'.csv'

#actual computation
for timeseries_index in range(8): #over time-series (previously manually inserted index of time-series to be used from timeseries variable)
    for blocknum in range(2,4): #over blocks (this means blocknum = 2 (1st round of the loop; blocknum = 3 (2nd round of the loop)))
        #loading data
        if blocknum == 2:
            #BLOCK 2 (2-8 Hz - OLD)
            mat = scipy.io.loadmat('/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/ToPython_2_8Hz/OLD.mat')
            timeseries = mat["lko"] #extracting timeseries
            blk = 'Old'
        elif blocknum == 3:
            #BLOCK 3 (2-8 Hz - NEW)
            mat = scipy.io.loadmat('/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/ToPython_2_8Hz/NEW.mat')
            timeseries = mat["lkn"] #extracting timeseries
            blk = 'New'
        
        if timeseries_index == 0 and blocknum == 2: #timeseries of parcel 1 for block OLD is non-fittable (which is fine since it is an occipital large brain parcel that has nearly nothing to do with the experimental task)
            lrow = '{} \n'
            lrow = lrow.format(timeseries_index + 1)
            csvfl.write(lrow)
        else:
                
            #fitting
            jes = list(range(time_point_start,time_point_end))
            data = timeseries[timeseries_index][jes]
            time = jes
            # time = list(range(1,jes[-1])) #array with time-points
            
            #defining function (summation of wavelets) for fitting
            def wavel(t, a, x0, sigma, c, phi, a2, x02, sigma2, c2, phi2, a3, x03, sigma3, c3, phi3, a4, x04, sigma4, c4, phi4, a5, x05, sigma5, c5, phi5):                 # define the function
                # return a*numpy.exp(-(t-x0)**2/(2*sigma**2))*numpy.cos(c*t+phi) + ho #normal gaussian
                #log gaussian (skewed)
                return (a*numpy.exp(-(numpy.log(t)-x0)**2/(2*sigma**2))*numpy.cos(c*t+phi) + a2*numpy.exp(-(numpy.log(t)-x02)**2/(2*sigma2**2))*numpy.cos(c2*t+phi2)
                        + a3*numpy.exp(-(numpy.log(t)-x03)**2/(2*sigma3**2))*numpy.cos(c3*t+phi3) + a4*numpy.exp(-(numpy.log(t)-x04)**2/(2*sigma4**2))*numpy.cos(c4*t+phi4)
                        + a5*numpy.exp(-(numpy.log(t)-x05)**2/(2*sigma5**2))*numpy.cos(c5*t+phi5))
             
            #BLOCK 2 - TIMESERIES 8,7,6,5,4,3,2, 1(more or less), NOT WITH 0 (very reasonable since it's an occipital big parcel)
            #BLOCK 3 - TIMESERIES 8,7,6,5,4,3,2, 1(more or less), 0(more or less)
            popt,pcov = curve_fit(wavel,time,data,p0=popo)
            # popt,pcov = curve_fit(wavel,time,data,p0=[-5,-1,0.05,-0.75,1.7, 12,3,0.7,0.1,3, -16,3,0.2,0.2,-0.4, -2.2,4,0.5,0.3,-7, 5,5,0.09,0.1,-7])
            
            # popt,pcov = curve_fit(wavel,time,data,p0=[0,0,3,0.19,0, 0,0,3,0.19,0, 0,0,1,0.19,0, 0,0,1,0.19,0, 0,0,1,0.19,0])
            a,x0,sigma,c,phi,a2,x02,sigma2,c2,phi2,a3,x03,sigma3,c3,phi3,a4,x04,sigma4,c4,phi4,a5,x05,sigma5,c5,phi5 = popt            
            #getting values predicted by the model        
            y_line = wavel(numpy.array(jes2),a,x0,sigma,c,phi,a2,x02,sigma2,c2,phi2,a3,x03,sigma3,c3,phi3,a4,x04,sigma4,c4,phi4,a5,x05,sigma5,c5,phi5) #getting predicted values
            data = timeseries[timeseries_index][jes2] #getting data of full time-series (for later plotting purposes)
    
            #getting R-squared value for model evaluation
            #residual sum of squares
            residuals = data - y_line 
            ss_res = numpy.sum(residuals**2)
            #total sum of squares
            ss_tot = numpy.sum((data-numpy.mean(data))**2)
            #R-squared value
            r_squared = 1 - (ss_res / ss_tot)
            
            #storing information in csv file
            lrow = '{};{};{};{};{};{};{};{};{};{};{};{};{};{};{};{};{};{};{};{};{};{};{};{};{};{};{};{} \n'
            lrow = lrow.format(timeseries_index + 1,blk,r_squared,a,x0,sigma,c,phi,a2,x02,sigma2,c2,phi2,a3,x03,sigma3,c3,phi3,a4,x04,sigma4,c4,phi4,a5,x05,sigma5,c5,phi5)
            csvfl.write(lrow)
            
            #plotting original data, subset of the data used for fitting function (if that is not the same as the full data), fitted function
            x_line = numpy.arange(min(timet),max(timet),1)
            # plt.plot(timet,timeseries[timeseries_index][0:315],color='orange') #plotting full timeseries
            plt.plot(timet,y_line,'--',color='red') #plotting fitted function (only within selected time-points)
            plt.plot(timet,data,color='blue') #plotting timeseries (only selected time-points)
            plt.ylim(-23,18) #y-axis limits
            plt.grid(b=True, which='both',linewidth=0.5) #grid
            plt.minorticks_on() #necessary to get minor grid
            #save figure
            pathfig = '/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/Figures/Figure_S7/Timeseries/N100_' + blk + '_parcel_' + str(timeseries_index + 1)
            plt.savefig(pathfig  + '.eps',format='eps')
            plt.savefig(pathfig  + '.png',format='png')
            plt.show()
            #saving timeseries for plotting them in Matlab (to get uniform graphics)
            savetxt(dir2 + 'fastfreq_predicted_data_' + str(timeseries_index + 1) + 'block_' + str(blocknum) + '.txt', y_line, delimiter = ',') #predicted data (y_line)
            savetxt(dir2 + 'fastfreq_real_data_' + str(timeseries_index + 1) + 'block_' + str(blocknum) + '.txt', data, delimiter = ',') #real data (data)

            
#closing csv file
csvfl.close()

#%%

#popo was calculated from a previous properly fitted function (similar to this one)
numpy.save('/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/ToPython_2_8Hz/popo.npy',popo)

#%%
