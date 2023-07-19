#/usr/bin/env ParselTongue

'''
ParselTongue runfile for LBA 22GHz data reduction

'''

from AIPS import AIPS
from definitions_parseltongue import *
import os, argparse

parser = argparse.ArgumentParser()
parser.add_argument("epoch",
                    help="Experiment code e.g. s006i",type=str)
parser.add_argument("aipsID",
                    help="AIPS ID, e.g 666",type=int)
parser.add_argument("-l","--log",
                    help="Logfile",type=str,default='ParselTongue.log')
parser.add_argument("-L","--linedisk",
                    help="Disk to save/find line data",type=int,default=1)
parser.add_argument("-r"."--refant",
                    help="Referance antenna",
                    type=int,default=0)
parser.add_argument("--load",
                    help="Load in FITS data",
                    action='store_true',default=False)
parser.add_argument("--prepr",
                    help="Calibrate line and cont data up \
                    up to manual phase cal",
                    action='store_true',default=False)
parser.add_argument("--mpcal",
                    help="Run manual phase cal",
                    action='store_true',default=False)
parser.add_argument("--cvel",
                    help="Run cvel",
                    action='store_true',default=False)
#parser.add_argument("--fring",
#                    help="Fringe fit on reference",
#                    action='store_true',default=False)
#parser.add_argument("--split",
#                    help="Split calibrators \
#                    and/or target channel",
#                    action='store_true',default=False)
#parser.add_argument("--image",
#                    help="Image calibrators \
#                    and/or target channel",
#                    action='store_true',default=False)

n       = 2   # number of uv-data files
defdisk = 1   # default disk
antname = 'LBA'

###############################################################
[filename, outname, outclass]    = [range(n),range(n),range(n)]
[nfiles, ncount, doconcat]       = [range(n),range(n),range(n)]
[outdisk, flagfile, antabfile]   = [range(n),range(n),range(n)]
for i in range(n):
    [flagfile[i], antabfile[i], outdisk[i]] = ['v490b.uvflag','',defdisk]
    [nfiles[i], ncount[i], doconcat[i]]     = [0,1,-1]
###############################################################

args        = parser.parse_args()
AIPS.userid = aipsID
logfile     = args.log

eop_path    = './eop/'
file_path   = '/path/to/fits/files/'

filename[0] = 'V490B.CONT.FITS'
outname [0] = 'V490B_CONT'
outclass[0] = 'UVDATA'
outdisk [0] = 1

filename[1] = 'V490B.LINE.FITS'
outname [1] = 'V490B_MASER'
outclass[1] = 'UVDATA'
outdisk [1] = args.linedisk

pos_shift   = {'N180' : [ -0.57, -0.53]}

cvelsource = ['IRAS00430-73','N66_SC12','IRAS01126-73','IRASF04521-6','NGC1984'
              ,'N113_MC24','N180','N157A_MC74','N160A_MC76','N214B'
              ,'HII_1107','HII_1186','IRAS05202_66','N105A_MC23','N159']

mp_source = ['0601-706']
mp_timera = [0,0,0,0,0,0,0,0]

cont   = 0
line   = 0
dpfour = 0

# only need to do this the first time
tasav_flag      = 1 

mprint('######################',logfile)
mprint(get_time(),logfile)
mprint('######################',logfile)

if args.load==true:
    for i in range(n):
        loadindx(file_path,filename[i],outname[i],outclass[i],outdisk[i],
                 nfiles[i],ncount[i],doconcat[i],antname,logfile)

data = range(n)
mprint('##################################',logfile)
for i in range(n):
    data[i]=AIPSUVData(outname[i], outclass[i], int(outdisk[i]), int(1))
    if data[i].exists():
        data[i].clrstat()
mprint('##################################',logfile)

for i in range(n):
    pr_data = data[i]

    if args.prepr==True:
        # restoring fg and su tables
        if tasav_flag==1: 
            runtasav(pr_data, i, logfile) 
        else:
            restore_fg(pr_data, logfile)
            restore_su(pr_data, logfile)
        # flaggin'
        runuvflg(pr_data,flagfile[i],logfile)
        # delete old sn and cl tables
        check_sncl(pr_data, 0, 1,logfile)
        # TEC CL1 -> CL2
        runTECOR(pr_data,year,doy,num_days,2,TECU_model)
        # EOPS CL2 -> CL3
        runeops(pr_data, eop_path)
        # PANG CL3 -> CL4
        runpang2(pr_data)
        
        # position shift CL4 -> CL4
        for source in pos_shift:
            [ra, dec] = [pos_shift[source][0],pos_shift[source][1]]
            if not source in pr_data.sources:
                continue
            if ra!=0 or dec!=0:
                if source == '':
                    source=findcal(pr_data, '')
                mprint('###############################'+
                       '###########################',logfile)
                mprint('Shift '+source+' by '+str(ra)+
                ' arcsec in RA and '+str(dec)+
                ' arcsec in DEC', logfile)
                mprint('###############################'+
                       '###########################',logfile)
                shift_pos(pr_data, source, ra, dec, 4, 4)
        
        # amplitude calibration:
        # make SN1 with accor
        runaccor(pr_data)
        # SN1 + CL4 = CL5
        runclcal(pr_data, 1, 4, 5, '',  1, args.refant)
        # CL5 -> CL6
        runapcal_lba(pr_data, 1, 5, 6)
    
    if args.mpcal==True:
        # make sure only up to SN1 and SN6 exist
        check_sncl(pr_data, 1, 6,logfile)
        (so,ti)=man_pcal(pr_data, args.refant, mp_source, mp_timera,debug,logfile, dpfour)
        if n==0:
            so_ti=[so,ti]
        if n==1:
            if so_ti[0]==so and so_ti[1]==ti:
                mprint('#############################################',logfile)
                mprint('### Both manual phasecal scans identical. ###',logfile)
                mprint('#############################################',logfile)
            else:
                mprint('#############################################',logfile)
                mprint('### Manual phasecal scans different.      ###',logfile)
                mprint('### Select one manually or close enough   ###',logfile)
                mprint('### '+data[1].name+' '+str(so_ti[0])+' '+str(so_ti[1])+' ###',logfile)
                mprint('### '+data[2].name+' '+str(so)+' '+str(ti)+' ###',logfile)
                mprint('#############################################',logfile)
                #sys.exit()
        # apply manual pcal, SN2 + CL6 -> CL7
        runclcal(pr_data, 2, 6, 7, '', 1, args.refant)        

if args.cvel==True:
    line_data = data[line]
    cont_data = data[cont]
    check_sncl(line_data, 2, 7,logfile)
    check_sncl(cont_data, 2, 7,logfile)
    runcvel(line_data,cvelsource,vel,inter_flag, doband, bpver)

