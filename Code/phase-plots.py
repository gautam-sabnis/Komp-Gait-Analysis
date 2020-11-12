import os 
import sys 
os.chdir("/Users/sabnig/Documents/Projects/Komp/Temp") 
import csv
import h5py
import math
import matplotlib
import matplotlib.pyplot as plt 
import pylab
import numpy as np
import pandas as pd
import seaborn as sns
import urllib.parse as urlparse
import gaitinference as ginf
import gaitplt as gplt

NUM_INTERP_FRAMES = 60
#NUM_INTERP_FRAMES = 2

cmap_tups = [
    (ginf.LEFT_FRONT_PAW_INDEX, 'Greens'),
    (ginf.RIGHT_FRONT_PAW_INDEX, 'Greens'),
    (ginf.LEFT_REAR_PAW_INDEX, 'Oranges'),
    (ginf.RIGHT_REAR_PAW_INDEX, 'Oranges'),
    (ginf.TIP_TAIL_INDEX, 'Blues'),
    (ginf.MID_TAIL_INDEX, 'Reds'),
    (ginf.BASE_TAIL_INDEX, 'Purples'),
    (ginf.NOSE_INDEX, 'Reds'),
    (ginf.BASE_NECK_INDEX, 'Blues'),
    (ginf.CENTER_SPINE_INDEX, 'Reds'), #'Greys'),
]
cmap_dict = dict(((str(k), v) for k, v in cmap_tups))

lines = [
    [ginf.TIP_TAIL_INDEX, ginf.MID_TAIL_INDEX, ginf.BASE_TAIL_INDEX,
     ginf.CENTER_SPINE_INDEX, ginf.BASE_NECK_INDEX, ginf.NOSE_INDEX],
    [ginf.LEFT_FRONT_PAW_INDEX, ginf.CENTER_SPINE_INDEX, ginf.RIGHT_FRONT_PAW_INDEX],
    [ginf.LEFT_REAR_PAW_INDEX, ginf.BASE_TAIL_INDEX, ginf.RIGHT_REAR_PAW_INDEX],
]


gait_h5 = h5py.File('../Data/KOMP-curated-2020-08-20.h5','r+')
metadata = pd.read_excel('../Data/KOMP-LinkedData.xlsx')

metadata = metadata.loc[metadata['NetworkFilename'].isin(urlparse.unquote(k) for k in gait_h5.keys()), :]
metadata['median_body_length_cm'] = float('nan')
for grp_name, grp in gait_h5.items():
    metadata.loc[metadata['NetworkFilename'] == urlparse.unquote(grp_name), 'median_body_length_cm'] = grp.attrs['median_body_length_cm']

stride_resolution = 100
def stride_bin_to_rad(stride_bin):
    return math.radians(3.6 * stride_bin)

speed_bin_size = gait_h5.attrs['speed_bin_size']
speed_bin_start = gait_h5.attrs['speed_bin_start']
speed_bin_stop = gait_h5.attrs['speed_bin_stop']
angular_velocity_bin_size = gait_h5.attrs['angular_velocity_bin_size']
angular_velocity_bin_count = gait_h5.attrs['angular_velocity_bin_count']

center_av_bin = -angular_velocity_bin_size // 2

# we only want to look at the centered angular velocity bin here
# angular_velocity_bin_count = 1

bin_count_dict = dict()
hildebrand_dict = dict()
norm_stride_pts_dict = dict()

# we need to segregate the files by strain
metadata = metadata.rename(columns = {"OFA_Genotype":"Strain","Mouse.ID":"MouseID"})
metadata.Strain = metadata.Strain.fillna('C57BL/6NJ')
metadata.Strain = [re.sub("<.*>","",x) for x in metadata.Strain]
metadata.Strain = [re.sub(" ","",x) for x in metadata.Strain]
metadata_grps = metadata.groupby('Strain')
for grp_key, metadata_grp in metadata_grps:
    speed_av_bins = ginf.gen_speed_and_av_bins(
        speed_bin_size, speed_bin_start,
        speed_bin_stop, angular_velocity_bin_size, angular_velocity_bin_count)

    for speed_av_bin in speed_av_bins:
        curr_speed, curr_av = speed_av_bin
        group_dict_key = (curr_speed, curr_av, grp_key,MouseID)
        bin_str = ginf.speed_av_bin_tup_to_str(speed_av_bin)
        for net_filename, mouse_id in zip(metadata_grp['NetworkFilename'], metadata_grp['MouseID']):
            escaped_file_name = urlparse.quote(net_filename, safe='')
            bin_path = escaped_file_name + '/bins/' + bin_str
            if bin_path in gait_h5:
                bin_grp = gait_h5[bin_path]
                both_paw_hild = np.stack([bin_grp['left_rear_hildebrand'], bin_grp['right_rear_hildebrand']])
                if group_dict_key not in bin_count_dict:
                    bin_count_dict[group_dict_key] = 1
                    hildebrand_dict[group_dict_key] = both_paw_hild
                    norm_stride_pts_dict[group_dict_key] = list(bin_grp['normalized_stride_points'])
                else:
                    bin_count_dict[group_dict_key] += 1
                    hildebrand_dict[group_dict_key] += both_paw_hild
                    norm_stride_pts_dict[group_dict_key] += bin_grp['normalized_stride_points']

ordered_keys = sorted(bin_count_dict.keys())

if True:
    tip_tail_lat_df_dict = dict()
    nose_lat_df_dict = dict()
    base_tail_lat_df_dict = dict()

    for curr_key in ordered_keys:
        speed, ang_vel, strain, MouseID = curr_key
        
        if ang_vel != center_av_bin:
            continue
        
        if speed not in (10,15,20,25,30):
            continue
        
        if strain not in ('Pcdh9-/+','C57BL/6NJ'):
        	continue

        print('=================================')
        print('speed: {}; strain: {}; count: {}'.format(
            speed, strain, bin_count_dict[curr_key]))

        curr_strides = ginf.restore_stride_points_shape(norm_stride_pts_dict[curr_key])
        interp_strides = np.stack(ginf.interpolate_stride_points(s, NUM_INTERP_FRAMES) for s in curr_strides)
        tip_tail_ys = interp_strides[:, :, 11, 1]
        tip_tail_ys -= np.repeat(np.mean(tip_tail_ys, axis=1), tip_tail_ys.shape[1]).reshape(tip_tail_ys.shape)
        nose_ys = interp_strides[:, :, 0, 1]
        nose_ys -= np.repeat(np.mean(nose_ys, axis=1), nose_ys.shape[1]).reshape(nose_ys.shape)
        base_tail_ys = interp_strides[:, :, 9, 1]
        base_tail_ys -= np.repeat(np.mean(base_tail_ys, axis=1), base_tail_ys.shape[1]).reshape(base_tail_ys.shape)
        num_strides = tip_tail_ys.shape[0]
        time_point = np.tile(np.arange(NUM_INTERP_FRAMES), num_strides)
        
        curr_tip_tail_lat_df = pd.DataFrame({
            'Percent Stride': np.tile(100.0 * np.arange(NUM_INTERP_FRAMES) / NUM_INTERP_FRAMES, num_strides),
            'Displacement': tip_tail_ys.flatten(),
            'stride_num': np.repeat(np.arange(num_strides), NUM_INTERP_FRAMES),
            'Mouse Line': strain,
            'Speed': speed,
        })
        curr_nose_lat_df = pd.DataFrame({
            'Percent Stride': np.tile(100.0 * np.arange(NUM_INTERP_FRAMES) / NUM_INTERP_FRAMES, num_strides),
            'Displacement': nose_ys.flatten(),
            'stride_num': np.repeat(np.arange(num_strides), NUM_INTERP_FRAMES),
            'Mouse Line': strain,
            'Speed': speed,
        })
        curr_base_tail_lat_df = pd.DataFrame({
            'Percent Stride': np.tile(100.0 * np.arange(NUM_INTERP_FRAMES) / NUM_INTERP_FRAMES, num_strides),
            'Displacement': base_tail_ys.flatten(),
            'stride_num': np.repeat(np.arange(num_strides), NUM_INTERP_FRAMES),
            'Mouse Line': strain,
            'Speed': speed,
        })

        tip_tail_lat_df_dict[(strain, speed)] = curr_tip_tail_lat_df
        base_tail_lat_df_dict[(strain, speed)] = curr_base_tail_lat_df
        nose_lat_df_dict[(strain, speed)] = curr_nose_lat_df


tip_tail_df = pd.DataFrame()
for k in tip_tail_lat_df_dict.keys():
	tip_tail_df = tip_tail_df.append(pd.DataFrame(tip_tail_lat_df_dict[k]))

nose_df = pd.DataFrame()
for k in nose_lat_df_dict.keys():
	nose_df = nose_lat_df.append(pd.DataFrame(nose_lat_df_dict[k]))

base_tail_df = pd.DataFrame()
for k in base_tail_lat_df_dict.keys():
    base_tail_df = base_tail_lat_df.append(pd.DataFrame(base_tail_lat_df_dict[k]))

strain_to_ctrl = {'Steap2-/-':'C57BL/6NJ'}

for mut_strain, ctrl_strain in strain_to_ctrl.items():
    print('Plots for mutant:', mut_strain)

    strains = (mut_strain, ctrl_strain)

    # here speed is fixed
    speed = 30

    ctrl_tip_tail_lat_dfs = [v for k, v in tip_tail_lat_df_dict.items() if k[0] in strains and k[1] == speed]
    ax = gplt.plot_lateral_disp(
        pd.concat(ctrl_tip_tail_lat_dfs, ignore_index=True),
        'Tip of Tail Lateral Displacement at {} cm/sec'.format(speed))
    plt.show()
    file_name = 'tip_tail_lat_disp_{0}.pdf'.format(re.sub("./.","",mut_strain))
    file_name = os.path.join('../Temp5/vignettes/', file_name)
    ax.get_figure().savefig(file_name, bbox_inches='tight')

    ctrl_nose_lat_dfs = [v for k, v in nose_lat_df_dict.items() if k[0] in strains and k[1] == speed]
    ax = gplt.plot_lateral_disp(
        pd.concat(ctrl_nose_lat_dfs, ignore_index=True),
        'Nose Lateral Displacement at {} cm/sec'.format(speed))
    plt.show()
    file_name = 'nose_lat_disp_{0}.pdf'.format(re.sub("./.","",mut_strain))
    file_name = os.path.join('../Temp5/vignettes/', file_name)
    ax.get_figure().savefig(file_name, bbox_inches='tight')

    ctrl_base_tail_lat_dfs = [v for k, v in base_tail_lat_df_dict.items() if k[0] in strains and k[1] == speed]
    ctrl_base_tail_lat_dfs = pd.concat(ctrl_base_tail_lat_dfs, ignore_index=True)
    ctrl_base_tail_lat_dfs['Speed'] = ctrl_base_tail_lat_dfs['Speed'].map('{} cm/sec'.format)
    ax = gplt.plot_lateral_disp(
        ctrl_base_tail_lat_dfs,
        'Base of Tail Lateral Displacement')
    plt.show()
    file_name = 'base_tail_lat_disp_{}.pdf'.format(re.sub("./.","",mut_strain))
    file_name = os.path.join('../Temp5/vignettes/', file_name)
    ax.get_figure().savefig(file_name, bbox_inches='tight')

    # here speed is included as hue
    ctrl_tip_tail_lat_dfs = [v for k, v in tip_tail_lat_df_dict.items() if k[0] in strains and k[1] == speed]
    ctrl_tip_tail_lat_dfs = pd.concat(ctrl_tip_tail_lat_dfs, ignore_index=True)
    ctrl_tip_tail_lat_dfs['Speed'] = ctrl_tip_tail_lat_dfs['Speed'].map('{} cm/sec'.format)
    ax = gplt.plot_lateral_disp2(
        ctrl_tip_tail_lat_dfs,
        'Tip of Tail Lateral Displacement')
    plt.show()
    file_name = 'tip_tail_lat_disp_{0}.pdf'.format(re.sub("./.","",mut_strain))
    file_name = os.path.join('../Temp5/vignettes/', file_name)
    ax.get_figure().savefig(file_name, bbox_inches='tight')

    ctrl_base_tail_lat_dfs = [v for k, v in base_tail_lat_df_dict.items() if k[0] in strains]
    ctrl_base_tail_lat_dfs = pd.concat(ctrl_base_tail_lat_dfs, ignore_index=True)
    ctrl_base_tail_lat_dfs['Speed'] = ctrl_base_tail_lat_dfs['Speed'].map('{} cm/sec'.format)
    ax = gplt.plot_lateral_disp2(
        ctrl_base_tail_lat_dfs,
        'Base of Tail Lateral Displacement')
    plt.show()
    file_name = 'base_tail_lat_disp_{}.pdf'.format(re.sub("./.","",mut_strain))
    file_name = os.path.join('../Temp5/vignettes/', file_name)
    ax.get_figure().savefig(file_name, bbox_inches='tight')


    ctrl_nose_lat_dfs = [v for k, v in nose_lat_df_dict.items() if k[0] in strains]
    ctrl_nose_lat_dfs = pd.concat(ctrl_nose_lat_dfs, ignore_index=True)
    ctrl_nose_lat_dfs['Speed'] = ctrl_nose_lat_dfs['Speed'].map('{} cm/sec'.format)
    ax = gplt.plot_lateral_disp2(
        ctrl_nose_lat_dfs,
        'Nose Lateral Displacement')
    plt.show()
    file_name = 'nose_lat_disp_{}.pdf'.format(re.sub("./.","",mut_strain))
    file_name = os.path.join('../Temp4/vignettes/', file_name)
    ax.get_figure().savefig(file_name, bbox_inches='tight')


#Base Tail Phase
base_tail_df = pd.concat([base_tail_df_Stmn4,base_tail_df_Med10,base_tail_df_Zfp579,base_tail_df_Sema6a,
    base_tail_df_Cmtr2], join = 'inner')

#df = base_tail_df[base_tail_df.Speed==30]
#ggplot(df, aes(x = 'Percent Stride', y = 'Displacement', color = 'Mouse Line')) + stat_smooth()

#R
#data <- read.csv('base_tail_df.csv', header=TRUE, stringsAsFactors=FALSE)
#df <- data[(data$Speed == 30) & (data$Mouse.Line %in% c('C57BL/6NJ','Stmn4-/-','Med10-/+','Zfp579-/-',
#    'Sema6a-/+','Cmtr2-/+')),]
#colnames(df)[5] <- 'Strain'
#ggplot(df, aes(x = Percent.Stride, y = Displacement, color = Strain)) + stat_smooth() + 
#labs(x = 'Percent Stride', y = 'Displacement (Base Tail)') + theme_bw(base_size=22) + theme(legend.position='top')
#ggsave('../Temp5/base_tail_vignette.pdf', width=9, height=9)


#Tip Tail Phase
tip_tail_df = pd.concat([tip_tail_df_Steap2,tip_tail_df_Dyx1c1,tip_tail_df_Ptcd3,tip_tail_df_Zfp579,
    tip_tail_df_Rai14], join = 'inner')
tip_tail_df.to_csv('tip_tail_df.csv')

#R
#data <- read.csv('tip_tail_df.csv', header=TRUE, stringsAsFactors=FALSE)
#df <- data[(data$Speed == 15) & (data$Mouse.Line %in% c('C57BL/6NJ','Steap2-/-','Dyx1c1-/+','Ptcd3-/+',
#    'Zfp579-/-','Rai14-/+')),]
#colnames(df)[5] <- 'Strain'
#ggplot(df, aes(x = Percent.Stride, y = Displacement, color = Strain)) + stat_smooth() + 
#labs(x = 'Percent Stride', y = 'Displacement (Tip Tail)') + theme_bw(base_size=22) + theme(legend.position='top')
#ggsave('../Temp5/tip_tail_vignette.pdf', width=9, height=9)

#Nose Phase
nose_df = pd.concat([nose_df_Steap2,nose_df_Prpf4b,nose_df_Rai14,nose_df_Stmn4,nose_df_Kcnd3], join = 'inner')
nose_df.to_csv('nose_df.csv')
#data <- read.csv('nose_df.csv', header=TRUE, stringsAsFactors=FALSE)
#df <- data[(data$Speed == 30),]
#colnames(df)[5] <- 'Strain'
ggplot(df, aes(x = Percent.Stride, y = Displacement, color = Strain)) + stat_smooth() + 
labs(x = 'Percent Stride', y = 'Displacement (Nose)') + theme_bw(base_size=22) + theme(legend.position='top')
#ggsave('../Temp5/nose_vignette.pdf', width=9, height=9)

