import matplotlib.pyplot as plt, pickle as pkl
from matplotlib import cm
import numpy as np
from sys import argv
from scipy import stats

script, pkl_file,sasa = argv

colors = cm.Dark2.colors
megalist = pkl.load(open(pkl_file,'rb'))
ifg = pkl_file.split('_')[0]

[water_contacts_ifglevel , hb_contacts_vdmlevel , vdw_contacts_vdmlevel,
    polar_contacts_vdmlevel, hb_contacts_atomlevel, ca_hb_contacts_atomlevel, 
    vdw_contacts_atomlevel, polar_contacts_atomlevel ] = megalist 
labels = ['# water molecules iFGs contact', '# vdMs iFGs hbond with', '# vdMs iFGs have vdw contacts with', '# vdMs iFGs have polar contacts with', '# atoms on vdMs that hbond with iFG', '# Ca atoms on vdMs that hbond with iFG', '# atoms on vdMs that make vdw contacts with iFG', '# atoms on vdMs that make polar contacts with iFG']  

########## REMOVE ELEMENTS IN LIST THAT ARE ABOVE SASA CUTOFF ##################




f, axarr = plt.subplots(4)
def plot(subplot,ix, each, cr):
    print(min(each), max(each))
    axarr[subplot].hist(each, bins=int(max(each))+1,histtype='step',color=colors[cr],label=labels[cr],lw=4)
    axarr[subplot].legend()
# top first 2 subplots, for each ifg (vdmlevel)
for ix, each in enumerate([water_contacts_ifglevel, vdw_contacts_vdmlevel]):
    plot(0,ix,each,ix)
for ix, each in enumerate([hb_contacts_vdmlevel, polar_contacts_vdmlevel]):
    plot(1,ix,each,ix+2)
# bottom subplot, for each vdm (atomlevel)
for ix, each in enumerate([hb_contacts_atomlevel, ca_hb_contacts_atomlevel]):
    flat = [item for sublist in each for item in sublist]
    plot(2,ix,flat,ix+4)
for ix, each in enumerate([vdw_contacts_atomlevel, polar_contacts_atomlevel]):
    flat = [item for sublist in each for item in sublist]
    plot(3,ix,flat,ix+6)
plt.suptitle('%s contacts - raw counts at %s A^2 cutoff' % (ifg,sasa))
plt.show()
