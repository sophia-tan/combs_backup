# assumes that ifg_interaction_types.pkl is in working dir
import matplotlib.pyplot as plt, pickle as pkl, pandas as pd
from matplotlib import cm
import numpy as np
from sys import argv
from scipy import stats

script, sasacsvpath, sasa = argv

colors = cm.Dark2.colors

sasacsv = pd.read_csv(sasacsvpath,index_col='iFG_count')
try:
    ifg = sasacsvpath.split('/')[7].split('_')[0]
except:
    ifg = sasacsvpath.split('/')[4].split('_')[0]

contactsdf = pkl.load(open('../setup_pkl/noskip_%s_contacts.pkl'%ifg,'rb'))
megalist = pkl.load(open('%s_interaction_types.pkl'%ifg,'rb'))

labels = ['# water molecules iFGs contact', '# vdMs iFGs have vdw contacts with', '# vdMs iFGs hbond with','# vdMs iFGs have polar contacts with', '# atoms on vdMs that hbond with iFG', '# Ca atoms on vdMs that hbond with iFG', '# atoms on vdMs that make vdw contacts with iFG', '# atoms on vdMs that make polar contacts with iFG']  

########## ONLY GET ELEMENTS IN LIST THAT ARE BELOW SASA CUTOFF ##################
# get iFGs that are below cutoff 
# poss_ifgs has more ifgs than is in the contacts pkl file bc contacts pkl file 
# has repeats removed
poss_ifgs = sasacsv.index[sasacsv['iFG_sasa_CB_3A_probe'] < float(sasa)].tolist()
# for every ifg in contactsdf (that has repeats removed), if it's also in poss_ifgs, then keep 
# the index number
contactsdf_ifgs = list(set(contactsdf['iFG_count']))
contactsdf_ifgs = np.array(contactsdf_ifgs)

######## get indices for ifgs that have 4 water contacts
#######unfiltered = [i for i, x in enumerate(megalist[0]) if x==4] 
######## use those indices to get ifg numbers
#######print(contactsdf_ifgs[unfiltered]) 

indices = np.isin(contactsdf_ifgs, np.array(poss_ifgs))
indices = list(np.where(indices==True)[0])
newlist = []
for ls in megalist: 
    newlist.append([ls[i] for i in indices])
megalist = newlist
[water_contacts_ifglevel , hb_contacts_vdmlevel , vdw_contacts_vdmlevel,
    polar_contacts_vdmlevel, hb_contacts_atomlevel, ca_hb_contacts_atomlevel, 
    vdw_contacts_atomlevel, polar_contacts_atomlevel ] = megalist 

########### only keep ifgs that pass cutoff so indices match ifg values
##########contactsdf_ifgs = contactsdf_ifgs[indices] 
########### get indices for ifgs that have 4 water contacts
##########unfiltered = [i for i, x in enumerate(water_contacts_ifglevel) if x==4]
########### use those indices to get ifg numbers
##########print(contactsdf_ifgs[unfiltered]) 



f, axarr = plt.subplots(4)
def plot(subplot,ix, each, cr):
    hist = np.histogram(each, bins=int(max(each))+1)
    x = np.arange(len(hist[0]))
    x = [i+ix/5 for i in x] # offset so histograms don't overlap
    axarr[subplot].bar(x,hist[0],color=colors[cr],label=labels[cr],lw=4,alpha=0.4,width=1)
    axarr[subplot].legend()
# top first 2 subplots, for each ifg (vdmlevel)
for ix, each in enumerate([water_contacts_ifglevel, vdw_contacts_vdmlevel]):
    plot(0,ix,each,ix)
for ix, each in enumerate([hb_contacts_vdmlevel, polar_contacts_vdmlevel]):
    plot(1,ix,each,ix+2)
# bottom subplot, for each vdm (atomlevel)
for ix, each in enumerate([hb_contacts_atomlevel, ca_hb_contacts_atomlevel]):
    print(len(each))
    flat = [item for sublist in each for item in sublist]
    plot(2,ix,flat,ix+4)
for ix, each in enumerate([vdw_contacts_atomlevel, polar_contacts_atomlevel]):
    flat = [item for sublist in each for item in sublist]
    plot(3,ix,flat,ix+6)
plt.suptitle('%s contacts - raw counts at %s A^2 cutoff' % (ifg,sasa))
#plt.show()
