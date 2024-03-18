
###########################################################
# 1. Inspect continuum of SgrB1off TM2
###########################################################

myvis = '/volumes/SHAOKC/2018.1.01420.S/calibrated_TM2/uid___A002_Xd58951_Xe45b.ms.split.cal'

# This will generate a text file that summarized the content of the data
# The text file is very useful, all the information we need below is taken from it.
listobs(vis=myvis, listfile=myvis+'.listobs')

# Set spws to be used (this information can be found in the listobs file)
contspws = '23,25,27,29'

# Split out SgrB1off pointing: #2 (this information can be found in the listobs file)
contvis = 'SgrB1off_TM2.ms'
split(vis=myvis, spw=contspws,      
    field='2', outputvis=contvis,
    datacolumn='all')

# This will generate another text file for the 'SgrB1off.ms' data
# Read this text file carefully
listobs(vis=contvis, listfile=contvis+'.listobs')

# Use plotms to identify lines and continuum
plotms(vis=contvis, xaxis='channel', yaxis='amplitude',
    field='0',
    uvrange='<30klambda',
    ydatacolumn='data',
    avgtime='1e8', avgscan=True, avgbaseline=True,
    iteraxis='spw',
    customsymbol=True,symbolshape='circle',symbolsize=4)

# Keep notes of line channels
# These are the channel indices with line meission.
# When we make continuum image, we need to exclude these channels
# But I didn't check it very carefully, and may miss some line emission.
# Please double check this part and let me know if you have any questions.
# The format is:
# - the number before the colon is the spectral window id (this information can be found in the listobs file)
# - the numbers in pairs after the colon are the channel ranges. Each pair marks a range of channels with line mission
# - the pairs are separated by semicolons
# - the spectral windows are separated by commas
fc = '0:412~450;1390~1460;1495~1520;1790~1830;3470~3550,1:480~500;520~740;3440~3530;3780~3820,2:150~250;700~800;1020~1070;1500~1670;3350~3430;3440~3490;3790~3830,3:370~480;1975~2080;3070~3200'

# Make a backup before flagging
flagmanager(vis=contvis,mode='save',
    versionname='before_cont_flags')

# Flag all the channels that have line emission
flagdata(vis=contvis,mode='manual',spw=fc,flagbackup=False)

# Check whether flags are as expected
plotms(vis=contvis,yaxis='amp',xaxis='channel',uvrange='<30klambda',
    avgtime='1e8',avgscan=True,iteraxis='spw',field='0',
    customsymbol=True,symbolshape='circle',symbolsize=4) 


###########################################################
# 2. Inspect continuum of SgrB1off TM1
###########################################################

myvis = '/volumes/SHAOKC/2018.1.01420.S/calibrated_TM1/uid___A002_Xe07f3e_Xf5d9.ms'##jump

# This will generate a text file that summarized the content of the data
# The text file is very useful, all the information we need below is taken from it.
listobs(vis=myvis, listfile=myvis+'.listobs')##jump

# Set spws to be used (this information can be found in the listobs file)
contspws = '23,25,27,29'

# Split out SgrB1off pointing: #2 (this information can be found in the listobs file)
contvis = 'SgrB1off_TM1_EB1.ms'
split(vis=myvis, spw=contspws,      
field='3', outputvis=contvis,
datacolumn='all')

# This will generate another text file for the 'SgrB1off.ms' data
# Read this text file carefully
listobs(vis=contvis, listfile=contvis+'.listobs')

# Use plotms to identify lines and continuum
plotms(vis=contvis, xaxis='channel', yaxis='amplitude',
    field='0',
    uvrange='<20klambda',
    ydatacolumn='data',
    avgtime='1e8', avgscan=True,
    iteraxis='spw',
    customsymbol=True,symbolshape='circle',symbolsize=4)

# Keep notes of line channels
# These are the channel indices with line meission.
# When we make continuum image, we need to exclude these channels
# But I didn't check it very carefully, and may miss some line emission.
# Please double check this part and let me know if you have any questions.
# The format is:
# - the number before the colon is the spectral window id (this information can be found in the listobs file)
# - the numbers in pairs after the colon are the channel ranges. Each pair marks a range of channels with line mission
# - the pairs are separated by semicolons
# - the spectral windows are separated by commas
fc = '0:360~500;1330~1600,1:550~730;3430~3530,2:150~230;740~820;1520~1660,3:360~500;1970~2100'###12345
# Make a backup before flagging
flagmanager(vis=contvis,mode='save',
    versionname='before_cont_flags')

# Flag all the channels that have line emission
flagdata(vis=contvis,mode='manual',spw=fc,flagbackup=False)

# check that flags are as expected.
plotms(vis=contvis,yaxis='amp',xaxis='channel',uvrange='<30klambda',
    avgtime='1e8',avgscan=True,iteraxis='spw',field='0',
    customsymbol=True,symbolshape='circle',symbolsize=4) 

###############jump116~172
# The same procedures, for the second EB in this SB
myvis = '/volumes/SHAOKC/2018.1.01420.S/calibrated_TM1/uid___A002_Xee9500_X286.ms'

# This will generate a text file that summarized the content of the data
# The text file is very useful, all the information we need below is taken from it.
listobs(vis=myvis, listfile=myvis+'.listobs')

# Set spws to be used (this information can be found in the listobs file)
contspws = '25,27,29,31'

# Split out SgrB1off pointing: #3 (this information can be found in the listobs file)
contvis = 'SgrB1off_TM1_EB2.ms'
split(vis=myvis, spw=contspws,      
field='3', outputvis=contvis,
datacolumn='all')

# This will generate another text file for the 'SgrB1off.ms' data
# Read this text file carefully
listobs(vis=contvis, listfile=contvis+'.listobs')

# Use plotms to identify lines and continuum
plotms(vis=contvis, xaxis='channel', yaxis='amplitude',
    field='0',
    uvrange='<30klambda',
    ydatacolumn='data',
    avgtime='1e8', avgscan=True,
    iteraxis='spw',
    customsymbol=True,symbolshape='circle',symbolsize=4)

# Keep notes of line channels
# These are the channel indices with line meission.
# When we make continuum image, we need to exclude these channels
# But I didn't check it very carefully, and may miss some line emission.
# Please double check this part and let me know if you have any questions.
# The format is:
# - the number before the colon is the spectral window id (this information can be found in the listobs file)
# - the numbers in pairs after the colon are the channel ranges. Each pair marks a range of channels with line mission
# - the pairs are separated by semicolons
# - the spectral windows are separated by commas
# EB2 spw2/3 were shifted by ~60 MHz, no idea why
fc = '0:400~460;1350~1550;1780~1850;3430~3620, 1:520~730;3430~3600;3760~3820, 2:10~100;630~685;1320~1660;3250~3300;3320~3380, 3:460~620;2050~2180;3200~3290'

# Make a backup before flagging
flagmanager(vis=contvis,mode='save',
    versionname='before_cont_flags')

# Flag all the channels that have line emission
flagdata(vis=contvis,mode='manual',spw=fc,flagbackup=False)

# check that flags are as expected.
plotms(vis=contvis,yaxis='amp',xaxis='channel',uvrange='<30klambda',
    avgtime='1e8',avgscan=True,iteraxis='spw',field='0',
    customsymbol=True,symbolshape='circle',symbolsize=4) 

###############################################
# Combine the two configurations
concat(vis=['SgrB1off_TM2.ms', 'SgrB1off_TM1_EB1.ms', 'SgrB1off_TM1_EB2.ms'], concatvis='SgrB1off_TM12_cont.ms')

myvis = 'SgrB1off_TM12_cont.ms'
listobs(vis=myvis, listfile=myvis+'.listobs')

## Make a smaller uv dataset
mstransform(vis='SgrB1off_TM12_cont.ms', outputvis='SgrB1off_TM12_cont_rebin30.ms', keepflags=False, chanaverage=True, chanbin=30, datacolumn='all')

# Rebin30: size from 712 GB to 22 GB
myvis = 'SgrB1off_TM12_cont_rebin30.ms'
listobs(vis=myvis, listfile=myvis+'.listobs')

# Recalculate the statistical weights
statwt(datacolumn='corrected', vis=myvis)

###############################################
# Image Parameters
cell = '0.04arcsec'
weighting = 'briggs'
robust = 0.5
niter = 100000000
# rms level is about 10 uJy/beam, and we set the threshold to ~3*rms
threshold = '0.03mJy'
imsize = [2560,2560]
imname = 'SgrB1off_band3_cont_rebin30v2'

tclean(vis = myvis,
  imagename = imname,
  specmode = 'mfs',
  deconvolver = 'mtmfs',
  nterms = 2,
  phasecenter = 'J2000 17:46:46.596292 -28.32.00.53396',
  scales = [0,5,15,50,150],
  imsize = imsize, 
  cell= cell, 
  gridder = 'mosaic',
  mosweight = True,
  usepointing = False,
  weighting = weighting, 
  robust = robust,
  niter = niter, 
  threshold = threshold, 
  interactive = True,
  pbcor = False,
  pblimit = 0.2)

impbcor(imagename=imname+'.image.tt0',pbimage=imname+'.pb.tt0',outfile=imname+'.image.tt0.pbcor')

exportfits(imagename=imname+'.image.tt0', fitsimage=imname+'.image.tt0.fits', dropdeg=True)
exportfits(imagename=imname+'.image.tt0.pbcor', fitsimage=imname+'.image.tt0.pbcor.fits', dropdeg=True)
exportfits(imagename=imname+'.pb.tt0', fitsimage=imname+'.pb.tt0.fits', dropdeg=True)
