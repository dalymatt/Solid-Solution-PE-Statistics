# -*- coding: utf-8 -*-
#Matthew Daly, Advanced Materials and Microstructures Lab, University of Illinois Chicago
#July 2022

#This function reads in the raw data from the Femtotools instrument and then
#produces a dataframe with the stored data. A second dataframe is produced with
#mapping information for the indents (positions, contact point)
#The function should be able to handle
#multiple indents in the same file, even those with varying conditions. But
#the indents should be all of the same type (static or CSM)
#flag=1, static indents
#flag=2, CSM indents
def import_raw_data(filename, flag=2):  

    import pandas as pd
    import numpy as np
    count=0
    
    if flag==1:
        dim=17
    elif flag==2:
        dim=27    
                
    with open(filename) as f:
        lines=f.readlines()
        dat=np.zeros((np.size(lines),dim))
        ind=np.ones((np.size(lines)))
        for line in lines:        
            if line.strip() and line[0]!='/' and line[1]!='/':
                myarray=np.fromstring(line, dtype=float, sep='\t')
                dat[count,:]=myarray[0:dim]
            else:
                ind[count]=0
            count=count+1

    dum=np.where(ind==0)[0]
    dum2=np.diff(dum)
    dum3=np.where(dum2>1)[0]
    dum3=dum3+1

    dum3=np.hstack((np.zeros(1),dum3))
    dum=np.hstack((dum, np.shape(dat)[0]))
    dum3=np.hstack((dum3, np.shape(dum)[0]-1))
    dum3=dum3.astype(int)

    indent_set=np.zeros((np.shape(dat)[0],2))
    
    for i in np.arange(0,np.shape(dum3)[0]-1):
        indent_set[dum[dum3[i]]:dum[dum3[i+1]],1]=i

    ind=ind.astype(bool)
    dat2=dat[ind,:]
    indent_set2=indent_set[ind,:]
  
    dum4=np.diff(dat2[:,3])
    dum5=np.where(dum4>=1)[0]
    
    dum5=dum5+1
    dum5=np.hstack((0,dum5))
    dum5=np.hstack((dum5, np.shape(dat2)[0]))
    k=0
 
    dum4=dum4.astype(int)
  
    for i in np.arange(0,np.shape(dum5)[0]-1):
        indent_set2[dum5[i]:dum5[i+1],0]=k
        k=k+1

    dat2=np.hstack((indent_set2,dat2))
   
    if flag==1:
       df=pd.DataFrame(dat2, columns=['Indent Num.', 'Indent set', 'Index [#]', 'Phase [#]', 'Displacement [um]', 'Time [s]', 'Pos X [um]',\
                                  'Pos Y [um]', 'Pos Z [um]', 'Rot A [m°]', 'Rot B [m°]', 'Piezo X [um]',\
                                      'Force A [uN]', 'Force B [uN]', 'Gripper [V]', 'Voltage A [V]',\
                                      'Voltage B [V]', 'Temperature [°C]',\
                                    'Sample Displace [um]'])
    elif flag==2:
        df=pd.DataFrame(dat2, columns=['Indent Num.', 'Indent set', 'Index [#]', 'Phase [#]', 'Displacement [um]', 'Time [s]', 'Pos X [um]',\
                                   'Pos Y [um]', 'Pos Z [um]', 'Rot A [m°]', 'Rot B [m°]', 'Piezo X [um]',\
                                       'Force A [uN]', 'Force B [uN]', 'Gripper [V]', 'Voltage A [V]',\
                                       'Voltage B [V]', 'Temperature [°C]', 'Amplitude force [uN]', 'Amplitude pos [um]',\
                                           'Stiffness [N/m]', 'Phase shift [°]', 'Phase excitation [°]',\
                                               'Force A raw [uN]', 'Displacement raw [um]', 'HMax [um]',\
                                                   'HMax raw [um]', 'Real Force [uN]', 'Real Stiffness [N/m]'])
  
    origin=df[['Pos Y [um]','Pos Z [um]']].values
    origin=origin[0,:]
    
    indent_array_map=np.zeros((int(np.max(dat2[:,0])+1),6))
    temp=np.zeros((1,1))

    for i in np.arange(0, np.shape(indent_array_map)[0]):        
    
        df_temp=df[df['Indent Num.']==i]        
        indent_array_map[i,0:3]=df_temp.head(1)[['Indent Num.', 'Indent set', 'Index [#]']].values
    
        indent_array_map[i,3]=df_temp.head(1)['Pos Y [um]'].values-origin[0]
        
        indent_array_map[i,4]=df_temp.head(1)['Pos Z [um]'].values-origin[1]
        
        if flag==1:
            displacement_offset=df_temp.head(1)['Displacement [um]'].values-df_temp.head(1)['Sample Displace [um]'].values
            indent_array_map[i,5]=displacement_offset
            piezo_offset=df_temp.head(1)['Piezo X [um]'].values-df_temp.head(1)['Sample Displace [um]'].values
            temp2=df_temp['Piezo X [um]'].values-piezo_offset
            
        elif flag==2:
            displacement_offset=df_temp.head(1)['Displacement [um]'].values-df_temp.head(1)['HMax [um]'].values
            indent_array_map[i,5]=displacement_offset
            piezo_offset=df_temp.head(1)['Piezo X [um]'].values-df_temp.head(1)['HMax [um]'].values
            temp2=df_temp['Piezo X [um]'].values-piezo_offset+ df_temp['Amplitude pos [um]'].values
        
        temp2=np.array([temp2])
        temp2=temp2.T
        temp=np.vstack((temp,temp2))
    temp=temp[1:,:]
    
    df['Piezo Sample Displace [um]']=temp
    df2=pd.DataFrame(indent_array_map, columns=['Indent Num.', 'Indent set', 'Index [#]','Rel. Pos Y [um]', 'Rel. Pos Z [um]', 'Contact point [um]'])    
 
    return df, df2
#%%
#This function will individually plot each indent's load displacement curve (without compliance correction)
#for both static and CSM indents. The static load signal is taken in both cases
def all_plot_raw(indents):

    import matplotlib as mpl
    import matplotlib.pyplot as plt        
    import numpy as np
    
    new_rc_params = {'text.usetex': False, "svg.fonttype": 'none'}
    mpl.rcParams.update(new_rc_params)
    
    def_font= {'fontname':'Arial'} #Setting a font name
    
    for i in np.arange(0,indents['Indent Num.'].max()+1):
        dat=indents[(indents['Indent Num.']==i)]
    
        plt.subplots(1,1, figsize=(8,5))
        ax1=plt.subplot(1,1,1)
        plt.plot(dat["Piezo Sample Displace [um]"],dat["Force A [uN]"],'k-')
    
        plt.yticks(**def_font, fontsize=16) 
        plt.xticks(**def_font, fontsize=16)
        plt.ylabel('Force A (uN)', fontsize=18, **def_font, labelpad=0)
        plt.xlabel('Piezo Sample Displace(um)', fontsize=18, **def_font, labelpad=0)
        plt.xlim(left=0)
        
        for axis in ['left', 'bottom', ]: #increasing the thickness of the plot border
            ax1.spines[axis].set_linewidth(1.5)
        for axis in ['top', 'right']: #increasing the thickness of the plot border
            ax1.spines[axis].set_visible(False)
        ax1.tick_params(axis='both', direction='in') #Setting the ticks to the inside of the plot
    return   
#%%
#This function will plot the location of indents on a planar map. It can accept
#both static and CSM indents
def array_map_plot(indents,indent_array_map):
 
    import matplotlib as mpl
    import matplotlib.pyplot as plt    
    import numpy as np
    import itertools    

    new_rc_params = {'text.usetex': False, "svg.fonttype": 'none'}
    mpl.rcParams.update(new_rc_params)
    
    def_font= {'fontname':'Arial'} #Setting a font name
        
    plt.subplots(1,1, figsize=(8,5))
    ax1=plt.subplot(1,1,1)
    marker = itertools.cycle(('s','v','^','d',',', '+', '.', 'o', '*'))
    
    for i in np.arange(0,indent_array_map['Indent set'].max()+1):        
        dat=indent_array_map[indent_array_map['Indent set']==i]        
        plt.plot(dat['Rel. Pos Y [um]'],dat['Rel. Pos Z [um]'],marker=next(marker), linestyle='')
    
    plt.yticks(**def_font, fontsize=16) 
    plt.xticks(**def_font, fontsize=16)
    plt.ylabel('Position z (um)', fontsize=18, **def_font, labelpad=0)
    plt.xlabel('Position y (um)', fontsize=18, **def_font, labelpad=0)
    
    for axis in ['left', 'bottom', ]: #increasing the thickness of the plot border
        ax1.spines[axis].set_linewidth(1.5)
    for axis in ['top', 'right']: #increasing the thickness of the plot border
        ax1.spines[axis].set_visible(False)
    ax1.tick_params(axis='both', direction='in') #Setting the ticks to the inside of the plot
 
    return 
#%%
#This function will separately plot indent sets, which should help identify outliers
#This function makes use of plotly so that rogue indents can be more easily identified
#The function can accept both CSM and static indent data.
#But, it does not guarantee the indents of the same set were performed with the same parameters
#Rather it just assumes that indents in the same set were collected with the same parameters 
def sep_plot(indents):
    import plotly.io as pio
    pio.renderers.default='browser'
    #import plotly.graph_objects as go
    import plotly.express as px
    
    for i in np.arange(0,int(indents['Indent set'].max()+1)):
        df=indents[indents['Indent set']==i]
        df=df[df['Piezo Sample Displace [um]']>0]
        #fig = go.Figure(data=[go.Line(x=df['Sample Displace [um]'], y=df['Force A [uN]'])])
        fig = px.line(df, x="Piezo Sample Displace [um]", y="Force A [uN]", color='Index [#]',\
                      custom_data=['Indent Num.','Indent set', 'Index [#]'])
            
        fig.update_traces(
        hovertemplate='Disp. [um] =%{x}<br> Load [uN] =%{y}<br> Indent Num.=%{customdata[0]}<br>Indent set =%{customdata[1]} <br> Index [#] =%{customdata[2]} <extra></extra>')
        fig.show()
    return
#%%
#This function will extract the unloading portions of load displacement curves.
#It should only be used for static indents. It will also return the compliance corrected
#indentation depth measurements.
#pick = an integer of the indent currently being examined
#Sf = the frame stiffness
#This function returns a dataframe with the unloading data with the compliance corrected
#indentation depth appended as a new column entry
def unload_extract(indents, pick, Sf=float('inf')):
 
    temp=indents[indents['Indent Num.']==pick]
    temp=temp[['Indent Num.', 'Indent set', 'Index [#]', 'Phase [#]','Piezo Sample Displace [um]', 'Force A [uN]']] 
  
    temp=temp[temp['Phase [#]']==4]
    
    temp['h sample corr [um]']=temp['Piezo Sample Displace [um]']-temp['Force A [uN]']/Sf
    
    a=temp[['h sample corr [um]']].idxmax()
    a=a.values
    a=a[0]
    temp=temp.loc[a:,:]
    
    unload=temp.copy()
    unload=unload.drop('Phase [#]',axis=1)
    unload.rename(columns={"Piezo Sample Displace [um]":"h sample [um]", "Force A [uN]": "P unload [uN]"}, inplace=True)
    cols=['Indent Num.', 'Indent set', 'Index [#]', 'h sample [um]', 'h sample corr [um]', 'P unload [uN]']
    unload=unload[cols]
   
    return unload
#%%
#This function will extract the loading portions of load displacement curves.
#It should only be used for CSM indents. It will return the Effective Stiffness,
#which is the  contact stiffness without compliance correction. This function achieves
#this by removing the compliance correction from the Femtotools stiffness measurements.
#pick = an integer of the indent currently being examined
#flag = 0 no plot of the frame stiffness measurement, 1 to plot it
#This function returns a dataframe with the loading data with the effective stiffness 
#appended as a new column entry. 

#Currently, the dynamic model for the Femtotools frame ignores
#a parallel spring for the indenter, and also we have ignored the mass of the indenter in dynamic calculations.
#See Hay 2010  "Continuous stiffness measurement during instrumented indentation testing" for more info
def load_extract(indents, pick, flag=0):

    import numpy as np    
    from scipy.optimize import curve_fit
    import matplotlib as mpl
    import matplotlib.pyplot as plt  
    
    temp=indents[indents['Indent Num.']==pick]
    temp=temp[['Indent Num.', 'Indent set', 'Index [#]', 'Phase [#]', 'Displacement [um]', 'Piezo X [um]',\
               'Piezo Sample Displace [um]','Force A [uN]', 'Real Force [uN]',\
                   'Amplitude force [uN]', 'Amplitude pos [um]', 'Phase shift [°]']] 
  
    temp=temp[temp['Phase [#]']==2]
    
    dat=(temp['Piezo X [um]'].values-temp['Piezo X [um]'].head(1).values)-\
        (temp['Displacement [um]'].values-temp['Displacement [um]'].head(1).values)
    dat2=temp['Force A [uN]'].values
    
    def stiffness_fit(x,m,b):
        return m*x+b
    
    params = curve_fit(stiffness_fit, dat, dat2, bounds=(-np.inf,np.inf), p0=np.array([0,dat2[0]]))
    params=params[0]
    Sf_NMT=params[0]
      
    temp3=(temp['Amplitude force [uN]'].values/temp['Amplitude pos [um]'].values)\
        *np.cos(np.radians(temp['Phase shift [°]'].values))
    
    if np.max(dat)>1e-5:
        temp3=(Sf_NMT*temp3)/(Sf_NMT+temp3)
    
    temp['Effective Stiffness [uN/um]']=temp3
    
    load=temp.copy()
    load=load.drop('Phase [#]',axis=1)
    load.rename(columns={"Piezo Sample Displace [um]":"h sample [um]", "Real Force [uN]": "P load [uN]"}, inplace=True)
    cols=['Indent Num.', 'Indent set', 'Index [#]', 'h sample [um]', 'Force A [uN]', 'P load [uN]',\
          'Amplitude force [uN]', 'Amplitude pos [um]', 'Phase shift [°]', 'Effective Stiffness [uN/um]']
    load=load[cols]   
    
    if flag==1:
        new_rc_params = {'text.usetex': False, "svg.fonttype": 'none'}
        mpl.rcParams.update(new_rc_params)
        
        def_font= {'fontname':'Arial'} #Setting a font name
        
        xfit=np.linspace(np.min(dat),np.max(dat),1000)
        yfit=stiffness_fit(xfit, params[0],params[1])
        
        plt.subplots(1,1, figsize=(8,5))
        ax1=plt.subplot(1,1,1)    
            
        plt.plot(dat,dat2,'ro', linestyle='', markerfacecolor='None')
        plt.plot(xfit,yfit,'k--')
        
        plt.yticks(**def_font, fontsize=16) 
        plt.xticks(**def_font, fontsize=16)
        plt.ylabel('Force A (uN)', fontsize=18, **def_font, labelpad=0)
        plt.xlabel('Piezo-displacement (um)', fontsize=18, **def_font, labelpad=0)
        
        for axis in ['left', 'bottom', ]: #increasing the thickness of the plot border
            ax1.spines[axis].set_linewidth(1.5)
        for axis in ['top', 'right']: #increasing the thickness of the plot border
            ax1.spines[axis].set_visible(False)
        ax1.tick_params(axis='both', direction='in') #Setting the ticks to the inside of the plot

    return load
#%%
#This function will extract the approach portion of the load depth curve. Currently
#this is not used, but may be useful later on.
def approach_extract(indents, pick, Sf=float('inf'), hcut=-1.0):
 
    temp=indents[indents['Indent Num.']==pick]
    temp=temp[['Indent Num.', 'Indent set', 'Index [#]', 'Phase [#]','Piezo Sample Displace [um]', 'Force A [uN]', 'Real Force [uN]',\
               'Amplitude force [uN]', 'Amplitude pos [um]', 'Phase shift [°]']] 
  
    temp=temp[temp['Phase [#]']==1]
    
    temp['h sample corr [um]']=temp['Piezo Sample Displace [um]']-temp['Force A [uN]']/Sf
    
    temp=temp[temp['h sample corr [um]']<=hcut]
    approach=temp.copy()
    approach=approach.drop('Phase [#]',axis=1)
    approach.rename(columns={"Piezo Sample Displace [um]":"h sample [um]", "Real Force [uN]": "P load [uN]"}, inplace=True)
    cols=['Indent Num.', 'Indent set', 'Index [#]', 'h sample [um]', 'h sample corr [um]', 'P load [uN]',\
          'Amplitude force [uN]', 'Amplitude pos [um]', 'Phase shift [°]']
    approach=approach[cols]   

    return approach

#%%
#This function will calculate the contact depth, contact stiffness at max load during unloading
#It does this by fitting the unloading curve between 100% and 50% max force to a polynomial expression.
#It will return a dataframe with these parameters
#If flag=1, the unloading curve and fitted expression will be plotted
def unload_params_static(unload, flag=0):

    import numpy as np
    import pandas as pd
    from scipy.optimize import curve_fit
    import matplotlib as mpl
    import matplotlib.pyplot as plt    
    
    Pup=1.0
    Plow=0.5
     
    dat=unload.values
    
    hmax=dat[0,4]
    Pmax=dat[0,5]
    
    ind=dat[:,5]<=Pmax*Pup
    ind2=dat[:,5]>=Pmax*Plow    
    ind3=ind2*ind
    
    dat2=dat[ind3,:]
    #m=1.5
    def OP_F(x,hf,m):
        return Pmax*((x - hf)/(hmax - hf))**m 
    
    unload_params_static=np.zeros((1,10))

    unload_params=curve_fit(OP_F, dat2[:,4], dat2[:,5], bounds=([0,1],[hmax,2]), p0=np.array([hmax/2,1.5])) 
    unload_params=unload_params[0]
    
    unload_params_static[0,0:5]=dat[0,[0,1,2,4,5]]
    
    unload_params_static[0,5]=unload_params[0]
    unload_params_static[0,6]=unload_params[1]
    S=(unload_params[1]*Pmax)/(hmax-unload_params[0])
    unload_params_static[0,7]=S
    epsilon=0.75
    unload_params_static[0,8]=hmax-(epsilon*(Pmax/S))
    residuals = dat2[:,5]- OP_F(dat2[:,4], *unload_params)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((dat2[:,5]-np.mean(dat2[:,5]))**2)
    r_squared = 1 - (ss_res / ss_tot)
    unload_params_static[0,9]=r_squared
 
    if flag==1:
        xfit=np.linspace(unload_params[0],hmax,100)
        yfit=Pmax*((xfit - unload_params[0])/(hmax - unload_params[0]))**unload_params[1]
            
        new_rc_params = {'text.usetex': False, "svg.fonttype": 'none'}
        mpl.rcParams.update(new_rc_params)
        
        def_font= {'fontname':'Arial'} #Setting a font name
        
        ind=100
        
        plt.subplots(1,1, figsize=(8,5))
        ax1=plt.subplot(1,1,1)    
            
        plt.plot(dat[::ind,4],dat[::ind,5],'ro', linestyle='', markerfacecolor='None')
        plt.plot(xfit,yfit,'k-')
        
        plt.yticks(**def_font, fontsize=16) 
        plt.xticks(**def_font, fontsize=16)
        plt.ylabel('Load (uN)', fontsize=18, **def_font, labelpad=0)
        plt.xlabel('Depth (um)', fontsize=18, **def_font, labelpad=0)
        plt.xlim(left=0)
        for axis in ['left', 'bottom', ]: #increasing the thickness of the plot border
            ax1.spines[axis].set_linewidth(1.5)
        for axis in ['top', 'right']: #increasing the thickness of the plot border
            ax1.spines[axis].set_visible(False)
        ax1.tick_params(axis='both', direction='in') #Setting the ticks to the inside of the plot
        
    unload_params_static=pd.DataFrame(unload_params_static, columns=['Indent Num.', 'Indent set', 'Index [#]','hmax [um]', 'Pmax [uN]', 'hf [um]',\
                                                       'm', 'Contact Stiffness [uN/um]', 'hc [um]', 'R^2'])

    return unload_params_static
#%%
#This function is used for CSM indents. It will return the compliance corrected indentation depth
#along with the compliance corrected contact stiffness. The contact stiffness is calculated from the
#effective stiffness using the expressions in Hay 2010 "Continuous stiffness measurement during instrumented indentation testing",
#With the same assumptions in the dynamic model (neglecting indenter mass and assuming no indenter spring in parallel to the contact+frame spring)
#This function will return a dataframe with the corrected indentation depth and contact stiffness, among other parameters.
#if flag is set to 1, the contact stiffness vs indentation depth will be plotted.
def stiffness_load_CSM_corr(load, Sf=float('inf'), flag=0):
 
    import numpy as np    
    import matplotlib as mpl
    import matplotlib.pyplot as plt    
    
    epsilon=0.75
    
    load['h sample corr [um]']=load['h sample [um]'].values-load['Force A [uN]'].values/Sf
    
    temp=load['Effective Stiffness [uN/um]'].values
    
    if Sf==np.inf:
        load['Contact Stiffness [uN/um]']=temp
    else :
        load['Contact Stiffness [uN/um]']=(Sf*temp)/(Sf-temp)
    
    load['hc [um]']=load['h sample corr [um]']-((load['P load [uN]'].values/load['Contact Stiffness [uN/um]'].values)*epsilon)
    
    if flag==1:
        new_rc_params = {'text.usetex': False, "svg.fonttype": 'none'}
        mpl.rcParams.update(new_rc_params)
        
        def_font= {'fontname':'Arial'} #Setting a font name
        
        plt.subplots(1,1, figsize=(8,5))
        ax1=plt.subplot(1,1,1)    
            
        plt.plot(load['h sample corr [um]'],load['Contact Stiffness corr [uN/um]'],'ro', linestyle='', markerfacecolor='None')
        
        plt.yticks(**def_font, fontsize=16) 
        plt.xticks(**def_font, fontsize=16)
        plt.ylabel('Contact Stiffness (uN/um)', fontsize=18, **def_font, labelpad=0)
        plt.xlabel('Distance into surface (um)', fontsize=18, **def_font, labelpad=0)
        
        for axis in ['left', 'bottom', ]: #increasing the thickness of the plot border
            ax1.spines[axis].set_linewidth(1.5)
        for axis in ['top', 'right']: #increasing the thickness of the plot border
            ax1.spines[axis].set_visible(False)
        ax1.tick_params(axis='both', direction='in') #Setting the ticks to the inside of the plot
    
    cols=['Indent Num.', 'Indent set', 'Index [#]', 'h sample [um]', 'h sample corr [um]', 'hc [um]', 'P load [uN]',\
          'Amplitude force [uN]', 'Amplitude pos [um]', 'Phase shift [°]', 'Effective Stiffness [uN/um]',\
              'Contact Stiffness [uN/um]']
    load=load[cols]   

    return load
#%%
#This function calculates the frame stiffness for static indents. It does this
#by assuming an array of stiffnesses and testing each one against the unloading data
#the stiffness that leads to the flattest P/S^2 curve will be selected
#flag is set to 1, two plots are provided. The first shows the P/S^2 vs indentation depth
#for a range of frame stiffnesses. The second plot shows the best fit
#The frame stiffness calibration data is output as a dataframe, with the
#best fit data output as a second dataframe
#P_S2_const is the P/S^2 constant of the material
#hcut is the indentation depth below which the P/S^2 data is ignored in fitting.
def Frame_stiffness_static(indents, flag=0,P_S2_const=0.0015,hcut=0.10):

    import numpy as np
    import pandas as pd
    from scipy.optimize import curve_fit
    import matplotlib as mpl
    import matplotlib.pyplot as plt  
    import itertools   
    
    indents_index=np.unique(indents['Indent Num.'].values)
    indents_index=indents_index.astype(int)       
    
    #An approximate range for the frame stiffness should be selected
    #This range should be estimated based on experience with different microforce sensors
    Sf_min=9800
    Sf_max=10500
    Sf_step=20
    # Sf_min=150000
    # Sf_max=180000
    # Sf_step=1000
    
    Sf_trial=np.arange(Sf_min,Sf_max,Sf_step)
    
    #Sf_trial=np.hstack((Sf_trial, np.inf))
    
    plot_ind=np.round(np.linspace(0,np.shape(Sf_trial)[0]-1,9).astype(int))
    
    def P_S2(x,slope):
        return slope*x+P_S2_const

    a=0

    for j in np.arange(0,np.shape(Sf_trial)[0]):

        print('Now at calibration step ' + str(j+1) + ' of ' + str(np.shape(Sf_trial)[0]))    
         
        indents_index_calib=indents_index
        for i in indents_index:
            unload=unload_extract(indents,i, Sf=Sf_trial[j])
            if unload.iloc[0,4] < hcut:
                ind=np.where(indents_index_calib==i)[0]
                indents_index_calib=np.delete(indents_index_calib,ind)
        if np.shape(indents_index_calib)[0]>0:
            k=0
            for i in indents_index_calib:
                unload=unload_extract(indents,i, Sf=Sf_trial[j])
                unload_params=unload_params_static(unload)
                if k==0:
                    unload_params_store=unload_params
                else:
                    unload_params_store=pd.concat([unload_params_store,unload_params])
                k=k+1
            Ps=unload_params_store['Pmax [uN]'].values/(unload_params_store['Contact Stiffness [uN/um]'].values\
                                                        *unload_params_store['Contact Stiffness [uN/um]'].values)
            Ps=Ps*1000
            unload_params_store['P/S^2 [GPa^-1]']=Ps.tolist()
        else:
            unload_params_store=pd.DataFrame({'Indent Num.':[np.nan],'Indent set': [np.nan], 'Index [#]': [np.nan],\
                                              'hmax [um]': [np.nan], 'Pmax [uN]': [np.nan], 'hf [um]': [np.nan],\
                                                  'm': [np.nan], 'Contact Stiffness [uN/um]': [np.nan], 'hc [um]': [np.nan], 'R^2': [np.nan],\
                                                      'P/S^2 [GPa^-1]':[np.nan]})
  
        unload_params_store = unload_params_store.reset_index(drop=True)
        
        Ps=unload_params_store['P/S^2 [GPa^-1]'].values
        
        if np.all(~np.isinf(Ps)) and np.all(~np.isnan(Ps)):
            params = curve_fit(P_S2, unload_params_store['hmax [um]'].values, Ps, bounds=(-1,1), p0=np.array([0]))
            params=params[0]            
            ss_res = Ps-P_S2(unload_params_store['hmax [um]'].values,*params)
            ss_res = np.sum(ss_res**2)
             
        else:
            params=np.nan            
            ss_res=np.nan            
                
        Calib_PS=pd.DataFrame(np.zeros((1,5)), columns=['Number of indents',  'Sf', 'Unload parameters',\
                                                             'Slope of P/S^2 [(GPa*um)^-1]', 'Sum of residual'])
        
        Calib_PS=Calib_PS.astype('object')
        Calib_PS.iloc[0,0]=np.shape(indents_index_calib)[0]
        Calib_PS.iloc[0,1]=Sf_trial[j]
               
        Calib_PS.iloc[0,2]=unload_params_store
        
        if ~np.isnan(params):
            Calib_PS.iloc[0,3]=params[0]
        else:
            Calib_PS.iloc[0,3]=np.nan
        Calib_PS.iloc[0,4]=ss_res
       
        if a==0:
            Calib_PS_store=Calib_PS
        else:
            Calib_PS_store=pd.concat([Calib_PS_store,Calib_PS])
            
        Calib_PS_store = Calib_PS_store.reset_index(drop=True)
        a=a+1
 
    b=np.nanargmin(np.abs((Calib_PS_store.iloc[:,3].values)*1000))
    
    Calib_use=Calib_PS_store.loc[[b]]
    Calib_use=Calib_use.reset_index(drop=True)
  
    if flag==1:
          new_rc_params = {'text.usetex': False, "svg.fonttype": 'none'}
          mpl.rcParams.update(new_rc_params)
          
          def_font= {'fontname':'Arial'} #Setting a font name
          
          marker = itertools.cycle(('s','v','^','d','>','+','<','p', 'o', '*'))      
          
          plt.subplots(1,1, figsize=(8,5))
          ax1=plt.subplot(1,1,1)    
          for i in np.arange(0,np.shape(plot_ind)[0]):        
              dat=Calib_PS_store.iloc[plot_ind[i],2].values        
              plt.plot(dat[:,3],(1000*dat[:,4])/(dat[:,7]*dat[:,7]),marker=next(marker), linestyle='')      
          
          plt.yticks(**def_font, fontsize=16) 
          plt.xticks(**def_font, fontsize=16)
          plt.ylabel('P/S^2 (GPa^-1)', fontsize=18, **def_font, labelpad=0)
          plt.xlabel('Corrected depth (um)', fontsize=18, **def_font, labelpad=0)
          plt.xlim(left=0)
          for axis in ['left', 'bottom', ]: #increasing the thickness of the plot border
              ax1.spines[axis].set_linewidth(1.5)
          for axis in ['top', 'right']: #increasing the thickness of the plot border
              ax1.spines[axis].set_visible(False)
          ax1.tick_params(axis='both', direction='in') #Setting the ticks to the inside of the plot
          
          temp=Calib_PS_store.iloc[-1,2].values
          
          xfit=np.linspace(0,np.max(temp[:,3]),100)
          yfit=xfit*Calib_use.iloc[0,3]+P_S2_const
          
          yfit2=P_S2_const*np.ones((np.shape(xfit)[0]))
          
          dat=Calib_use.iloc[0,2].values
            
          plt.subplots(1,1, figsize=(8,5))
          ax1=plt.subplot(1,1,1)    
        
          plt.plot(dat[:,3],(1000*dat[:,4])/(dat[:,7]*dat[:,7]),'ro')
          plt.plot(xfit,yfit,'k-')
          plt.plot(xfit,yfit2,'b--')
           
          plt.yticks(**def_font, fontsize=16) 
          plt.xticks(**def_font, fontsize=16)
          plt.ylabel('P/S^2 (GPa^-1)', fontsize=18, **def_font, labelpad=0)
          plt.xlabel('Corrected depth (um)', fontsize=18, **def_font, labelpad=0)
          plt.xlim(left=0)
          plt.ylim([0, 0.01])
          for axis in ['left', 'bottom', ]: #increasing the thickness of the plot border
              ax1.spines[axis].set_linewidth(1.5)
          for axis in ['top', 'right']: #increasing the thickness of the plot border
              ax1.spines[axis].set_visible(False)
          ax1.tick_params(axis='both', direction='in') #Setting the ticks to the inside of the plot      
   
    return Calib_use, Calib_PS_store
#%%
#This function calculates the frame stiffness for CSM indents. It does this
#by assuming an array of stiffnesses and testing each one against the loading data
#the stiffness that leads to the flattest P/S^2 curve will be selected
#flag is set to 1, two plots are provided. The first shows the P/S^2 vs indentation depth
#for a range of frame stiffnesses. The second plot shows the best fit
#The frame stiffness calibration data is output as a dataframe, with the
#best fit data output as a second dataframe
#P_S2_const is the P/S^2 constant of the material
#hcut is the indentation depth below which the P/S^2 data is ignored in fitting.
def Frame_stiffness_CSM(indents, flag=0,P_S2_const=0.0015,hcut=0.05):
    
    import numpy as np
    import pandas as pd
    from scipy.optimize import curve_fit
    import matplotlib as mpl
    import matplotlib.pyplot as plt  
    import itertools   
    
    indents_index=np.unique(indents['Indent Num.'].values)
    indents_index=indents_index.astype(int)       
    
    #An approximate range for the frame stiffness should be selected
    #This range should be estimated based on experience with different microforce sensors
    Sf_min=9000
    Sf_max=11000
    Sf_step=10
    # Sf_min=1e5
    # Sf_max=2e5
    # Sf_step=1000
    
    Sf_trial=np.arange(Sf_min,Sf_max,Sf_step)
    
    Sf_trial=np.hstack((Sf_trial, np.inf))
    
    plot_ind=np.round(np.linspace(0,np.shape(Sf_trial)[0]-1,9).astype(int))
    
    def P_S2(x,slope):
        return slope*x+P_S2_const

    a=0

    for j in np.arange(0,np.shape(Sf_trial)[0]):
  
        print('Now at calibration step ' + str(j+1) + ' of ' + str(np.shape(Sf_trial)[0]))    

        k=0    
        for i in indents_index:
            load=load_extract(indents,i, flag=0)
            load_corr=stiffness_load_CSM_corr(load, Sf=Sf_trial[j])
            if k==0:
                load_corr_store=load_corr
            else:
                load_corr_store=pd.concat([load_corr_store,load_corr])
            k=k+1
        
        Ps=load_corr_store['P load [uN]'].values/(load_corr_store['Contact Stiffness [uN/um]'].values\
                                                  *load_corr_store['Contact Stiffness [uN/um]'].values)
        Ps=Ps*1000
        load_corr_store['P/S^2 [GPa^-1]']=Ps
        
        load_corr_store = load_corr_store.reset_index(drop=True)
        
        temp=load_corr_store['h sample corr [um]']>=hcut
        
        params = curve_fit(P_S2, load_corr_store[temp]['h sample corr [um]'].values, Ps[temp], bounds=(-1,1), p0=np.array([0]))
        params=params[0]            
        ss_res = Ps[temp]-P_S2(load_corr_store[temp]['h sample corr [um]'].values,*params)
        ss_res = np.sum(ss_res**2)
             
        Calib_PS=pd.DataFrame(np.zeros((1,5)), columns=['Number of indents',  'Sf', 'Load parameters',\
                                                             'Slope of P/S^2 [(GPa*um)^-1]', 'Sum of residual'])
        
        Calib_PS=Calib_PS.astype('object')
        Calib_PS.iloc[0,0]=np.shape(indents_index)[0]
        Calib_PS.iloc[0,1]=Sf_trial[j]
               
        Calib_PS.iloc[0,2]=load_corr_store
                
        Calib_PS.iloc[0,3]=params[0]
        
        Calib_PS.iloc[0,4]=ss_res
       
        if a==0:
            Calib_PS_store=Calib_PS
        else:
            Calib_PS_store=pd.concat([Calib_PS_store,Calib_PS])
            
        Calib_PS_store = Calib_PS_store.reset_index(drop=True)
        
        a=a+1
 
    b=np.nanargmin(np.abs((Calib_PS_store.iloc[:,3].values)*1000))
    
    Calib_use=Calib_PS_store.loc[[b]]
    Calib_use=Calib_use.reset_index(drop=True)
    
    if flag==1:
          new_rc_params = {'text.usetex': False, "svg.fonttype": 'none'}
          mpl.rcParams.update(new_rc_params)
          
          def_font= {'fontname':'Arial'} #Setting a font name
          
          marker = itertools.cycle(('s','v','^','d','>','+','<','p', 'o', '*'))      
          
          plt.subplots(1,1, figsize=(8,5))
          ax1=plt.subplot(1,1,1)    
          for i in np.arange(0,np.shape(plot_ind)[0]):        
              dat=Calib_PS_store.iloc[plot_ind[i],2].values        
              plt.plot(dat[:,4],dat[:,12],marker=next(marker), linestyle='')      
          
          plt.yticks(**def_font, fontsize=16) 
          plt.xticks(**def_font, fontsize=16)
          plt.ylabel('P/S^2 (GPa^-1)', fontsize=18, **def_font, labelpad=0)
          plt.xlabel('Corrected depth (um)', fontsize=18, **def_font, labelpad=0)
          plt.xlim(left=0)
          for axis in ['left', 'bottom', ]: #increasing the thickness of the plot border
              ax1.spines[axis].set_linewidth(1.5)
          for axis in ['top', 'right']: #increasing the thickness of the plot border
              ax1.spines[axis].set_visible(False)
          ax1.tick_params(axis='both', direction='in') #Setting the ticks to the inside of the plot
          
          temp=Calib_PS_store.iloc[-1,2].values
          
          xfit=np.linspace(0,np.max(temp[:,4]),100)
          yfit=xfit*Calib_use.iloc[0,3]+P_S2_const
          
          yfit2=P_S2_const*np.ones((np.shape(xfit)[0]))
          
          dat=Calib_use.iloc[0,2].values
            
          plt.subplots(1,1, figsize=(8,5))
          ax1=plt.subplot(1,1,1)    
        
          plt.plot(dat[:,4],dat[:,12],'ro')
          plt.plot(xfit,yfit,'k-')
          plt.plot(xfit,yfit2,'b--')
           
          plt.yticks(**def_font, fontsize=16) 
          plt.xticks(**def_font, fontsize=16)
          plt.ylabel('P/S^2 (GPa^-1)', fontsize=18, **def_font, labelpad=0)
          plt.xlabel('Corrected depth (um)', fontsize=18, **def_font, labelpad=0)
          plt.xlim(left=0)
          plt.ylim([0, 0.004])
          for axis in ['left', 'bottom', ]: #increasing the thickness of the plot border
              ax1.spines[axis].set_linewidth(1.5)
          for axis in ['top', 'right']: #increasing the thickness of the plot border
              ax1.spines[axis].set_visible(False)
          ax1.tick_params(axis='both', direction='in') #Setting the ticks to the inside of the plot    
    
    return Calib_use, Calib_PS_store
#%%
#This function is used to remove outlier indents from the datasets. It returns
#a dataframe with the truncated raw data
def exclude_indents(indents,exclude=[]):
    drops=indents_data[indents_data['Indent Num.'].isin(exclude)].index
    indents_trunc=indents_data.drop(drops)
    return indents_trunc
#%%
#This function will return the coefficients for the tip function of a Berkovich tip,
#following the Oliver Pharr model. #The calibration material is assumed to be fused silica.
#The reduced modulus is can be changed (Er_fused)
#Beta is the tip correction factor used in the Oliver Pharr model. Typical range is between
#1-1.05 (see Oliver 2004, "Measurement of hardness and elastic modulus by instrumented 
#indentation: Advances in understanding and refinements to methodology" )
#if flag is set to 1, the tip function, hardness, reduced modulus, and P/S^2 curves
#are plotted against contact depth. 
#This function will work with both CSM and static indentatino data

def Tip_function_fit(Calib_use,Er_fused=69.6, beta=1.034, P_S2_const=0.0015, flag=0):

    import numpy as np
    from scipy.optimize import curve_fit
    import pandas as pd
    import matplotlib as mpl
    import matplotlib.pyplot as plt  
    
    # E_fused=72
    # E_ind=1141
    # nu_ind=0.07
    # nu_fused=0.17
    # Er_fused=((1-(nu_ind**2))/E_ind)+((1-(nu_fused**2))/E_fused)
    # Er_fused=1/Er_fused
    
    # beta=1.034 
    
    Hardness_fused=0.0015*(4*beta*beta)*Er_fused*Er_fused/np.pi
    
    calibrated=Calib_use.iloc[0,2]
    
    temp=(np.pi/4)*(((calibrated['Contact Stiffness [uN/um]'].values)/(beta*Er_fused))**2)/(1e6)
    calibrated['A [um^2]']=temp.tolist()
    
    def A_OP(x,C0,C1,C2,C3,C4,C5,C6,C7,C8):
        return (C0*x*x)+(C1*x)+(C2*(x**(1/2)))+(C3*(x**(1/4)))+(C4*(x**(1/8)))+\
            (C5*(x**(1/16)))+(C6*(x**(1/32)))+(C7*(x**(1/64)))+(C8*(x**(1/128)))
    
    p0s=np.array([24.5,0,0,0,0,0,0,0,0])
    boundsLow=np.array([20,0,0,0,0,0,0,0,0])
    boundsUp=np.array([30,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf])
    
    params_A_OP=curve_fit(A_OP, calibrated['hc [um]'].values,\
                               calibrated['A [um^2]'].values, bounds=(boundsLow,boundsUp), p0=p0s)
    params_A_OP=params_A_OP[0]
    
    ss_res = calibrated['A [um^2]'].values-A_OP(calibrated['hc [um]'].values,*params_A_OP)
    ss_res = np.sum(ss_res**2)
    
    ss_tot = np.sum(((calibrated['A [um^2]'].values-np.mean(calibrated['A [um^2]'].values))**2))
    Tip_R_squared = 1 - (ss_res / ss_tot)
    
    Er=np.sqrt((np.pi/(4*beta*beta))*(calibrated['Contact Stiffness [uN/um]'].values*calibrated['Contact Stiffness [uN/um]'].values)*\
        (1/A_OP(calibrated['hc [um]'].values,*params_A_OP)))/1000
    
    H=Er*Er*calibrated['P/S^2 [GPa^-1]'].values*((4*beta*beta)/np.pi)
    
    calibrated['E reduced [GPa]']=Er
    calibrated['Hardness [GPa]']=H
    
    Tip_out = np.hstack((params_A_OP,Tip_R_squared))
    Tip_out=np.array([Tip_out])
    Tip_function=pd.DataFrame(Tip_out, columns=['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6',\
                                                         'C7', 'C8','R^2'])        
    if flag==1:
          new_rc_params = {'text.usetex': False, "svg.fonttype": 'none'}
          mpl.rcParams.update(new_rc_params)
          
          def_font= {'fontname':'Arial'} #Setting a font name
          
          plt.subplots(1,1, figsize=(8,5))
          ax1=plt.subplot(1,1,1)    
                    
          plt.plot(calibrated['hc [um]'],calibrated['A [um^2]'], 'ro')      
    
          xfit=np.linspace(0.01,calibrated['hc [um]'].max(),100)
          yfit=A_OP(xfit,*params_A_OP)            
                  
              
          plt.plot(xfit,yfit,'k--')      
          
          plt.yticks(**def_font, fontsize=16) 
          plt.xticks(**def_font, fontsize=16)
          plt.ylabel('A (um^2)', fontsize=18, **def_font, labelpad=0)
          plt.xlabel('hc (um)', fontsize=18, **def_font, labelpad=0)
          plt.xlim(left=0)
          plt.ylim(bottom=0)
          for axis in ['left', 'bottom', ]: #increasing the thickness of the plot border
              ax1.spines[axis].set_linewidth(1.5)
          for axis in ['top', 'right']: #increasing the thickness of the plot border
              ax1.spines[axis].set_visible(False)
          ax1.tick_params(axis='both', direction='in') #Setting the ticks to the inside of the plot
          
          
          def consts_plot(x,Er_fused):
              return Er_fused+0*x
          
          yfit_Er=consts_plot(xfit,Er_fused)
          
          plt.subplots(1,1, figsize=(8,5))
          ax1=plt.subplot(1,1,1)    
                    
          plt.plot(calibrated['hc [um]'],calibrated['E reduced [GPa]'], 'cv', markerfacecolor='none')      
          plt.plot(xfit,yfit_Er,'k--')      
          
          plt.yticks(**def_font, fontsize=16) 
          plt.xticks(**def_font, fontsize=16)
          plt.ylabel('Er (GPa)', fontsize=18, **def_font, labelpad=0)
          plt.xlabel('hc (um)', fontsize=18, **def_font, labelpad=0)
          plt.xlim(left=0)
          plt.ylim(bottom=0)
          plt.ylim(top=calibrated['E reduced [GPa]'].max()+30)
          for axis in ['left', 'bottom', ]: #increasing the thickness of the plot border
              ax1.spines[axis].set_linewidth(1.5)
          for axis in ['top', 'right']: #increasing the thickness of the plot border
              ax1.spines[axis].set_visible(False)
          ax1.tick_params(axis='both', direction='in') #Setting the ticks to the inside of the plot
          
          yfit_H=consts_plot(xfit,Hardness_fused)
          
          plt.subplots(1,1, figsize=(8,5))
          ax1=plt.subplot(1,1,1)    
                    
          plt.plot(calibrated['hc [um]'],calibrated['Hardness [GPa]'], 'gs', markerfacecolor='none')      
          plt.plot(xfit,yfit_H,'k--')      
          
          plt.yticks(**def_font, fontsize=16) 
          plt.xticks(**def_font, fontsize=16)
          plt.ylabel('Hardness (GPa)', fontsize=18, **def_font, labelpad=0)
          plt.xlabel('hc (um)', fontsize=18, **def_font, labelpad=0)
          plt.xlim(left=0)
          plt.ylim(bottom=0)
          plt.ylim(top=calibrated['Hardness [GPa]'].max()+3)
          for axis in ['left', 'bottom', ]: #increasing the thickness of the plot border
              ax1.spines[axis].set_linewidth(1.5)
          for axis in ['top', 'right']: #increasing the thickness of the plot border
              ax1.spines[axis].set_visible(False)
          ax1.tick_params(axis='both', direction='in') #Setting the ticks to the inside of the plot
          
          yfit_P_S2=consts_plot(xfit,P_S2_const)
          
          plt.subplots(1,1, figsize=(8,5))
          ax1=plt.subplot(1,1,1)    
                    
          plt.plot(calibrated['hc [um]'],calibrated['P/S^2 [GPa^-1]'], 'b^', markerfacecolor='none')      
          plt.plot(xfit,yfit_P_S2,'k--')      
          
          plt.yticks(**def_font, fontsize=16) 
          plt.xticks(**def_font, fontsize=16)
          plt.ylabel('P/S^2 (GPa^-1)', fontsize=18, **def_font, labelpad=0)
          plt.xlabel('hc (um)', fontsize=18, **def_font, labelpad=0)
          plt.xlim(left=0)
          plt.ylim(bottom=0)
          plt.ylim(top=calibrated['P/S^2 [GPa^-1]'].max()+0.002)
          for axis in ['left', 'bottom', ]: #increasing the thickness of the plot border
              ax1.spines[axis].set_linewidth(1.5)
          for axis in ['top', 'right']: #increasing the thickness of the plot border
              ax1.spines[axis].set_visible(False)
          ax1.tick_params(axis='both', direction='in') #Setting the ticks to the inside of the plot
          
    return calibrated, Tip_function


#%%
#This function wll analyze indentation data, assuming a tip function is known. If 
#the frame stiffness is not known, it will assume a value of infinity. 
#The function will output the projected contact area, the reduced modulus, the elastic modulus,
#and the hardness of the indented material
#This function accepts both CSM and static indentation data. The elastic properties of the indenter
#can be proivded along with the Poisson's ratio of the indented material.
#flag = 1 is for static indents
#flag = 2 is for CSM indents
def Indent_analyze(indents,flag=1,E_ind=1141,nu_ind=0.07,nu_mat=0.17, Sf_use=float('inf')):

    import numpy as np
    import pandas as pd

    E_ind=1141
    nu_ind=0.07
    nu_mat=0.17
    
    indents_index=np.unique(indents['Indent Num.'].values)
    indents_index=indents_index.astype(int)
    
    k=0
    for i in indents_index:
        if flag==1:
            dat=unload_extract(indents,i, Sf=Sf_use)
            dat_params=unload_params_static(dat)
        elif flag==2:
            dat=load_extract(indents,i)
            dat_params=stiffness_load_CSM_corr(dat, Sf=Sf_use)

        if k==0:
            indents_analyze=dat_params
        else:
            indents_analyze=pd.concat([indents_analyze,dat_params])
        k=k+1

    def Area_function(x,C0,C1,C2,C3,C4,C5,C6,C7,C8):
        return (C0*x*x)+(C1*x)+(C2*(x**(1/2)))+(C3*(x**(1/4)))+(C4*(x**(1/8)))+\
            (C5*(x**(1/16)))+(C6*(x**(1/32)))+(C7*(x**(1/64)))+(C8*(x**(1/128)))
    
    tip_function_coeffs=Tip_function.iloc[0,:-1].values
    
    indents_analyze['Area [um^2]']=Area_function(indents_analyze['hc [um]'].values,*tip_function_coeffs)
    
    Er=(indents_analyze['Contact Stiffness [uN/um]'].values)*(np.sqrt(np.pi)/2)*(1/np.sqrt(indents_analyze['Area [um^2]']))/1000
    indents_analyze['Er [GPa]']=Er
    
    E_mat=1/Er
    
    E_mat=E_mat-((1-nu_ind*nu_ind)/E_ind)
    
    E_mat=1/E_mat
    
    E_mat=E_mat*(1-nu_mat*nu_mat)
    
    indents_analyze['E [GPa]']=E_mat
    if flag==1:
        indents_analyze['Hardness [GPa]']=indents_analyze['Pmax [uN]']/indents_analyze['Area [um^2]']/1000
    elif flag==2:
        indents_analyze['Hardness [GPa]']=indents_analyze['P load [uN]']/indents_analyze['Area [um^2]']/1000
        
    return indents_analyze
      
#%%
import numpy as np

np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning) 

filename='CSM_Data_FusedSI.txt' #Change to name of raw text file from nanoindenter

#Read indentation data into one dataframe and provide characteristics of indentation array
#flag = 1, static indents
#flag = 2, CSM indents
indents_data,indent_array_map=import_raw_data(filename, flag=2) 

#all_plot_raw(indents_data) #Plot each indent individually

#array_map_plot(indents_data,indent_array_map) #Plots a map of the indents

#sep_plot(indents_data) #Plot indentation replications on same plot and organize replica indents into lists

#indents_data_trunc=exclude_indents(indents_data,exclude=[5]) Exclude specific indents (see sep_plot to choose)

#sep_plot(indents_data_trunc) #Plot indentation replications on same plot and organize replica indents into lists

#Stiffness calibration for static indents
#flag = 0 Do not produce plots
#flag = 1 Produce P/S^2 plots
#Calib_use, Calib_static_store=Frame_stiffness_static(indents_data, flag=1) 

#Stiffness calibration for CSM indents
#flag = 0 Do not produce plots
#flag = 1 Produce P/S^2 plots
Calib_use, Calib_CSM_store=Frame_stiffness_CSM(indents_data, flag=1)

#Tip function fitting
#flag = 0, Do not plot Area vs contact depth, H vs contact depth, Er vs contact depth plots, or P/S^2 vs contact depth plots
#flag = 1, Plot Area vs contact depth, H vs contact depth, Er vs contact depth plots, or P/S^2 vs contact depth plots
Indent_calibrated, Tip_function=Tip_function_fit(Calib_use,flag=1)

#Calculate indentation properties (Area, Er, E, H) for indents with a calibrated tip function
#flag = 1, Static indents
#flag = 2 CSM indents
#Sf = frame stiffness
Indents_analysis=Indent_analyze(indents_data, Sf_use=float(Calib_use['Sf'].values), flag=2)




