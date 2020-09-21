#!/usr/bin/env python
from beamformlib import *
import glob

bool_tab = True # 1 if coherent beamforming, 0 if incoherent beamforming
bool_aips = 1 # 1 if AIPS calibration, 0 if MIR calibration
bool_excludeak01 = 0 # 1 if exclude ak01, 0 if not
frb = '200430' #"191001" #vela191001   180924, 181112, 190102, 190608, 190611.2, 190711, 191001, 191228

# VCRAFT files
basedr = "/fred/oz002/users/hcho/craft/"
#files = glob.glob(basedr+"python/voltages/FRB"+frb+"/*/*/*.vcraft")
#files = glob.glob("/fred/oz002/users/sbhandari/rawdata_frb191001/vela/*/*/*.vcraft") # Vela directory
files = glob.glob("/fred/oz002/askap/craft/frb200430/frb/rawdata/SB13593/20200430154900/ICS/C000/ak**/beam4[89]/*.vcraft")
if bool_excludeak01: # exclude ak01
    new_files =[]
    for f in files:
        if "ak01" not in f: # or "ak33" not in f:
            new_files += [f]
    files = new_files
files_x,files_y = sort_by_pol(files)

i = 32768 #1 #2048*3
n = 1
outfile = "corr.fits"
an = None
if frb == "181112":
    # FRB181112 details
    if bool_aips: # AIPS calibration
        calcfile = basedr+"Calibration/aipscal/frb181112/c1_f0/craftfrb.im"
        #"Calibration/mircal/frb181112/SB7030_FRB_neweop.im"
        hwfile = basedr+"Calibration/mircal/frb181112/hwdelays_0407_SB7031_20181112222901_round-8.txt"
        parset = basedr+"Calibration/aipscal/frb181112/fcm.txt"
        mir = None
        aips = basedr+"Calibration/aipscal/frb181112/bandpass.bp.txt"
    else:
        calcfile = basedr+"Calibration/mircal/frb181112/SB7030_FRB_neweop.im"
        hwfile = basedr+"Calibration/mircal/frb181112/hwdelays_0407_SB7031_20181112222901_round-8.txt"
        #"Calibration/FRB181112/SB7031oldeop.hwdelays"
        parset = basedr+"Calibration/mircal/frb181112/fcm_from_odin_181112.txt"
        #"Calibration/FRB180924-calcfiles/fcm_release2_ak25mod.txt"
        mir = basedr+"Calibration/mircal/frb181110/20181112222901_call_beam00_i4096_f9.uvlin"
        #"Calibration/0407_allfreq/20180924212734_call_beam37_i1024_f9.uvaver" #None
        aips = None
    offset = 10003+32*56000+32*9500 #+32*30000

if frb == "180924":
    # FRB180924 details
    print("FRB180924")
    offset=1874193
    calcfile=basedr+"Calibration/aipscal/frb180924/c1_f0/craftfrb.im" # geometric delays
    parset=basedr+"Calibration/aipscal/frb180924/fcm.txt" # clock delays
    hwfile=basedr+"Calibration/SB6635_b18_neweop.hwdelays" # hardware delays
    mir=basedr+"Calibration/0407_allfreq/20180924212734_call_beam37_i1024_f9.uvaver" # MIR gain, bandpass
    aips=basedr+"Calibration/aipscal/frb180924/bandpass.bp.txt" # AIPS gain, bandpass
    tag="_aips"
    
if frb == "190102":
    # FRB190102 details
    print("FRB190102")
    offset=1880759 #1750775
    calcfile=basedr+"Calibration/aipscal/frb190102/craftfrb.im" # geometric delays
    parset=basedr+"Calibration/aipscal/frb190102/fcm.txt" # clock delays
    hwfile=None # hardware delays
    if not bool_tab:
        aips=None
        mir=None # MIR gain, bandpass
    else:
        aips=basedr+"Calibration/aipscal/frb190102/noxpol/bandpasses.bp.txt" # AIPS gain, bandpass
        mir=None
        tag="_aips"

if frb == "190608":
    # FRB190608 details
    print("FRB"+frb)    
    i = 16384 #temp
    n = 1
    offset=2018819 #1750775
    calcfile=basedr+"Calibration/aipscal/frb190608/craftfrb.im" # geometric delays
    parset=basedr+"Calibration/aipscal/frb190608/fcm.txt" # clock delays
    hwfile=None # hardware delays
    aips=basedr+"Calibration/aipscal/frb190608/noxpol/bandpasses.bp.txt" # AIPS gain, bandpass
    mir=None
    tag="_aips"

if frb == "190611.2":
    # FRB190611.2 details
    print("FRB"+frb)    
    i = 16384
    n = 1
    offset=1887164
    calcfile=basedr+"Calibration/aipscal/frb190611/craftfrb.im"#_200418.im" # geometric delays
    parset=basedr+"Calibration/aipscal/frb190611/fcm.txt" # clock delays
    hwfile=None # hardware delays
    if not bool_tab:
        aips=None
    else:
        aips=basedr+"Calibration/aipscal/frb190611/noxpol/bandpasses_noxpol_FRB190611.2.bp.txt" # AIPS gain, bandpass
        tag="_aips"
    mir=None

if frb == "190711":
    # FRB190711 details
    print("FRB"+frb)    
    i = 16384
    n = 1
    offset=302100 #242868
    calcfile="/fred/oz002/askap/craft/frb190711/FRB/rficor/data/c5_f0/craftfrb_8544.im" #basedr+"Calibration/aipscal/frb190714/craftfrb.im"#_200418.im" # geometric delays
    parset="/fred/oz002/askap/craft/frb190711/FRB/rficor/data/c5_f0/fcm.txt" #basedr+"Calibration/aipscal/frb190714/fcm.txt" # clock delays
    hwfile=None # hardware delays
    if not bool_tab:
        aips=None
    else:
        aips=None #basedr+"Calibration/aipscal/frb190611/noxpol/bandpasses_noxpol_FRB190611.2.bp.txt" # AIPS gain, bandpass
        tag="_aips"
    mir=None

if frb == "191001":
    # FRB191001 details
    print("FRB"+frb)    
    #i = 57408 #Full 3 sec of data
    #i = 32768 # 1.8sec of data
    i = 25926 # 1.4 sec of data	
    #i = 22222 # 1.2 sec of data
    n = 1
    #offset=2251851 # Skipping 1.9 sec of data. 
    #offset= 2464 #3344 #3344 #48 #3344 # 39984
    offset=2014832 # Skipping 1.7 sec of data
    calcfile="/fred/oz002/adeller/askap/frb191001/model-for-hyerin/craftfrb_182569.im"
    #"/fred/oz002/askap/craft/frb191001/FRB/data.gate/c1_f0/craftfrb_182569.im" # geometric delays
    parset="/fred/oz002/askap/craft/frb191001/FRB/data.gate/c1_f0/fcm.txt" # clock delays
    hwfile=None # hardware delays
    if not bool_tab:
        aips=None
    else:
        aips="/fred/oz002/askap/craft/frb191001/FRB/data.gate/bandpasses_noxpol.bp.txt" # AIPS gain, bandpass
        tag="_aips"
    mir=None

if frb == "vela191001":
    # Beamforming on vela
    print("FRB"+frb)
    #i = 57408 #Full 3 sec of data
    #i = 32768 # 1.8sec of data
    i = 57408
    n = 1
    offset=3372
    #offset= 2464 #3344 #3344 #48 #3344 # 39984
    calcfile="/fred/oz002/askap/craft/frb191001/vela/data.vela.gate/c1_f0/craftfrb_185647.im"
    #"/fred/oz002/askap/craft/frb191001/FRB/data.gate/c1_f0/craftfrb_182569.im" # geometric delays
    parset="/fred/oz002/askap/craft/frb191001/vela/data.vela.gate/c1_f0/fcm.txt" # clock delays
    hwfile=None # hardware delays
    if not bool_tab:
        aips=None
    else:
        aips="/fred/oz002/askap/craft/frb191001/vela/data.vela.gate/bandpasses_noxpol.bp.txt" # AIPS gain, bandpass
        tag="_aips"
    mir=None


if frb == "191228":
    # FRB191228 details
    print("FRB"+frb)    
    i = 16384
    n = 1
    offset=1576224 #3360
    calcfile="/fred/oz002/askap/craft/frb191228/frb/newreftime/data.1ms.cube.nrt/c6_f0/craftfrb_148602.im" # geometric delays
    parset="/fred/oz002/askap/craft/frb191228/frb/newreftime/data.1ms.cube.nrt/c6_f0/fcm.txt" # clock delays
    hwfile=None # hardware delays
    if not bool_tab:
        aips=None
    else:
        aips=basedr+"Calibration/aipscal/frb"+frb+"/noxpol/bandpasses_noxpol.bp.txt" # AIPS gain, bandpass
        tag="_aips"
    mir=None


if frb == '200430':
    offset=230
    DM=381.045
    f0=864.5
    #calcfile="/fred/oz002/users/dscott/Calibration/aipscal/frb200430/craftfrb_139162.im"
    calcfile="/fred/oz002/askap/craft/frb200430/frb/data.1ms.dmshift2/c1_f0/craftfrb_204531.im"
    #fcm="/fred/oz002/users/dscott/Calibration/aipscal/frb200430/fcm.txt"
    fcm="/fred/oz002/askap/craft/frb200430/frb/data.1ms.dmshift2/c1_f0/fcm.txt"
    hwfile=
    mir=
    aips="/fred/oz002/users/dscott/Calibration/aipscal/frb200430/bandpasses_noxpol_FRB200430.bp.txt"

    if [ "$pol" == "x" ]; then
        f_vcraft="/fred/oz002/askap/craft/frb200430/frb/rawdata/SB13593/20200430154900/ICS/C000/ak**/beam48/*.vcraft"
    elif [ "$pol" == "y" ]; then
        f_vcraft="/fred/oz002/askap/craft/frb200430/frb/rawdata/SB13593/20200430154900/ICS/C000/ak**/beam49/*.vcraft"
    else
        echo "ERROR: Must provide polarisation as x or y!"
    fi

    i=32768
    n=1
    n_ant=26



# 'x' polarization
values = ParseValues(i,n,calcfile,hwfile,parset,mir,aips,an,offset,files_x,outfile,bool_tab=bool_tab)

corr = get_corr(values)
if bool_tab:
    sum_tot = np.power(abs(corr.do_tab(values.an)),2)
else:
    sum_tot = abs(corr.do_tab(values.an))

# save x polarization
if bool_aips:
    tag = "_aips"
else:
    tag = "_mir"
if bool_excludeak01:
    tag += "_excludeak01"

if not bool_tab:
    bf_type = "inco"
    tag = ""
else:
    bf_type = "tab"
outfile = 'test_output/beamformed/'+bf_type+'_frb'+frb+'_n'+str(n)+'_i'+str(i)+tag+'_polx.npy'
print('X polarization intensity saved to: '+outfile)
np.save(outfile,sum_tot)

# 'y' polarization
#files = glob.glob(basedr+"python/voltages/FRB181112/**/beam01/*.vcraft")
#values = ParseValues(i,n,calcfile,hwfile,parset,mir,aips,an,offset,files,outfile)
if bool_aips == 0:
    mir = basedr+"Calibration/mircal/frb181112/20181112222901_call_beam01_i4096_f9.uvlin"
    values.change_mir(mir)
values.change_files(files_y)


corr = get_corr(values)

if bool_tab:
    sum_tot += np.power(abs(corr.do_tab(values.an)),2)
else:
    sum_tot += abs(corr.do_tab(values.an))

# save total intensity
outfile = 'test_output/beamformed/'+bf_type+'_frb'+frb+'_n'+str(n)+'_i'+str(i)+tag+'.npy'
print('saving to: '+outfile)
np.save(outfile,sum_tot)
