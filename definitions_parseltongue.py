import numpy as np

from AIPS import AIPS
from AIPSTask import AIPSTask, AIPSList
from AIPSData import AIPSUVData, AIPSImage
from Wizardry.AIPSData import AIPSUVData as WAIPSUVData

##############################################################################
# TEC stuff
##############################################################################

def get_TEC(year,doy,TECU_model):
    year=str(year)[2:4]
    if doy<10:
        doy='00'+str(doy)
    elif doy<100:
        doy='0'+str(doy)
    else:
        doy=str(doy)
    name=TECU_model+doy+'0.'+year+'i'
    if os.path.exists(name):
        print 'File already there.'
    else:
        #path='ftp://cddis.gsfc.nasa.gov/gps/products/ionex/20'+year+'/'+doy+'/'
        path='ftp://gdc.cddis.eosdis.nasa.gov/gnss/products/ionex/20'+year+'/'+doy+'/'
        #os.popen(r'wget -t 30 -O '+name+'.Z '+path+name+'.Z')
        os.popen(r'curl --insecure -O --ftp-ssl '+path+name+'.Z')
        os.popen(r'uncompress -f '+name+'.Z')

def runTECOR(indata,year,doy,num_days,gainuse,TECU_model):
    year=str(year)[2:4]
    if doy<10:
        doy='00'+str(doy)
    if doy<100:
        doy='0'+str(doy)
    else:
        doy=str(doy)
    name=TECU_model+doy+'0.'+year+'i'
    tecor = AIPSTask('TECOR')
    if os.path.exists(name):
        tecor.infile='PWD:'+name
    tecor.indata=indata
    tecor.nfiles=num_days
    tecor.gainuse = gainuse
    tecor.aparm[1:] = [1,0]
    tecor() 

##############################################################################
# EOP stuff
##############################################################################

def get_eop(eop_path):
    if os.path.exists(eop_path+'usno_finals.erp'):
        age = (time.time() - os.stat(eop_path+'usno_finals.erp')[8])/3600
        mprint('usno_finals.erp exists, not downloaded.',logfile)
    else:
        os.popen(r'curl --insecure -O --ftp-ssl ftp://gdc.cddis.eosdis.nasa.gov/vlbi/gsfc/ancillary/solve_apriori/usno_finals.erp')   #wget ftp://cddis.gsfc.nasa.gov/vlbi/gsfc/ancillary/solve_apriori/usno500_finals.erp http://gemini.gsfc.nasa.gov/solve_save ftp://ftp.lbo.us/pub/staff/wbrisken/EOP
        os.popen(r'curl --insecure -O --ftp-ssl ftp://gdc.cddis.eosdis.nasa.gov/vlbi/gsfc/ancillary/solve_apriori/usno500_finals.erp')
        os.popen(r'mv usno500_finals.erp '+eop_path+'usno_finals2.erp')
        os.popen(r'mv usno_finals.erp '+eop_path+'usno_finals.erp')

#
def runeops(indata, eop_path, cl_in=2, cl_out=3):
    eops        = AIPSTask('CLCOR')
    eops.indata = indata
    eops.gainver = cl_in
    eops.gainuse = cl_out
    eops.opcode  = 'EOPS'
    eops.infile  = eop_path+'usno_finals.erp'
    eops()

##############################################################################
# Pang
##############################################################################

def runpang2(indata):
    pang        = AIPSTask('CLCOR')
    pang.indata = indata
    pang.gainver = 3
    pang.gainuse = 4
    pang.opcode  = 'PANG'
    pang.clcorprm[1:] = [1,0]
    antennas          = []
    # the next bit is to deal with HB,KE,YG being linear
    for row in indata.table('AN', 0):
        if row['mntsta']==0:
            if not row['anname'].replace(' ','') in ('HB','KE'):
                antennas.append(row['nosta'])
    pang.antennas[1:] = antennas
    pang()

##############################################################################
# Flagging
##############################################################################

def runuvflg(indata,flagfile,logfile):
    if flagfile!='' and os.path.exists(flagfile):
        uvflg        = AIPSTask('UVFLG')
        uvflg.indata = indata
        uvflg.intext = flagfile
        uvflg.opcode = 'FLAG'
        uvflg.go()
    else:
        mprint('No UVFLG file applied.',logfile)

##############################################################################
# Loading
##############################################################################

def loadindx(filepath,filename,outname,outclass,
             outdisk,nfiles,ncount,doconcat,antname,logfile):
    if os.path.exists(filepath+filename):
        mprint('File exists!',logfile)
    else:
        raise RuntimeError('File does not exists!')

    fitld = AIPSTask('FITLD')
    fitld.datain   = filepath+filename
    fitld.outname  = outname
    fitld.outclass = outclass
    fitld.outseq   = 1
    fitld.outdisk  = int(outdisk)
    fitld.ncount   = ncount
    fitld.nfiles   = nfiles
    fitld.doconcat = doconcat
    fitld.clint    = 1./60.
    fitld.wtthresh = 0.45
    if aipsver!='31DEC09':
        if not antname=='LBA': 
            fitld.antname[1:] = [antname]

    data = AIPSUVData(fitld.outname, fitld.outclass,
                      int(fitld.outdisk), int(fitld.outseq))
    if data.exists():
        data.clrstat()
        data.zap()
        mprint('##############################',logfile)
        mprint('Data already there => deleted!',logfile)
        mprint('##############################',logfile)
    else:
        mprint('#########################',logfile)
        mprint('Data not there => read in',logfile)
        mprint('#########################',logfile)

    fitld.go()

    mprint('################################################',logfile)
    mprint(str(data)+' loaded!',logfile)
    mprint('################################################',logfile)

    if data.exists():
        data.zap_table('AIPS CL',1)
        runindxr(data)
        mprint('#################',logfile)
        mprint('Data new indexed!',logfile)
        mprint('#################',logfile)
    else:
        mprint('No!',logfile)
    return

##############################################################################
# Amplitude calibration
##############################################################################

 def runaccor(indata):
    accor           = AIPSTask('ACCOR')
    accor.indata    = indata
    accor.timer[1:] = [0]
    accor.solint    = 2
    accor()

# 22GHz 

def runapcal_lba(indata, snver, inver, outver):
    sqrtsefd        = {'AT':10.49, 'CD': 44.72, 'HH': 38.73, 
                      'HO': 44.72, 'MP': 30.00, 'PA':14.14}
    ant             = get_ant(indata)
    check_sncl(indata, snver, inver, logfile)
    for i in ant:
        clcor             = AIPSTask('CLCOR')
        clcor.opcode      = 'GAIN'
        clcor.indata      = indata
        clcor.gainver     = 0
        clcor.gainuse     = outver
        clcor.antenna[1:] = [i]
        clcor.clcorprm[1] = sqrtsefd[ant[i]]
        clcor.go()
    return

##############################################################################
# Velocity corrections
##############################################################################

def get_center_freq(indata):
    fq = indata.table('AIPS FQ',0)
    naxis = indata.header['naxis']
    if naxis[3]>1:
        fq_span=fq[0]['if_freq'][indata.header.naxis[3]-1]-fq[0]['if_freq'][0]
        frq=(indata.header.crval[2]+0.5*fq_span)
    else:
        frq=indata.header.crval[2]
    return frq

def get_line_name(indata):
    freq = get_center_freq(indata)/1e9
    if freq>12 and freq<13:
        restfreq = [1.2178E+10,597000]
        linename='CH3OH_12GHz'
        print 'Assuming 12.2 GHz methanol maser.'
    elif freq>22 and freq<23:
        restfreq = [2.2235E+10,80000]
        linename='H2O'
        print 'Assuming 22.2 GHz water maser.'
    elif freq>6 and freq<7:
        restfreq = [6.668e+09, 519200]
        linename='CH3OH_6.7GHz'
        print 'Assuming 6.7 GHz methanol maser.'
    elif freq>43.1 and freq<43.2:
        restfreq = [43.122e+09,80000]
        linename='SiO_43.1GHz'
        print 'Assuming 43.122 GHz SiO maser.'
    else:
        print 'Unknown maser line.'
        exit()
    return (linename, restfreq)

def findcal(indata, calsource):
    if calsource == '':
        n = 0
        for source in indata.sources:
            if source[0]=='G':
                calsource=source
                n=n+1
        if n>1:
            print 'More than one Maser source! Using '+calsource
    return calsource

def run_setjy(indata, source, flux):
    setjy = AIPSTask('SETJY')
    setjy.indata = indata
    setjy.source[1:] = [source]
    setjy.zerosp[1:] = flux
    setjy.optype = ''
    setjy.bif    = 0
    setjy.eif    = 0
    setjy.optype = 'VCAL'
    print setjy.optype
    setjy()

def findcvelsource(indata, cvelsource):
    if cvelsource == '' or cvelsource == ['']:
        cvelsource = []
        n = 0
        for source in indata.sources:
            if source[0]=='G' and source[1:3]!='RB':
                cvelsource.append(source)
                n=n+1
    return cvelsource

def runcvel_lba(indata, cvelsource, inter_flag, doband, bpver):
    #uvsort(indata)
    print 'Running CVEL.'
    if inter_flag==1:
        cvelsource=check_calsource(indata, cvelsource)

    naxis = indata.header['naxis']
    crpix = indata.header['crpix']

    (linename,restfreq) = get_line_name(indata)

    setjy = AIPSTask('SETJY')
    setjy.indata = indata
    setjy.source[1:] = cvelsource
    setjy.restfreq[1:] = restfreq
    setjy.optype = 'VCAL'
    setjy.veltyp = 'LSR'
    setjy.veldef = 'RADIO'
    channum = indata.header['naxis'][2]
    setjy.aparm[1:] = [channum/2.+crpix[2],0]

    for i in range(naxis[3]):
        #setjy.sysvel     = vel[i]
        setjy.bif = i+1
        setjy.eif = i+1
        setjy()

    cvel = AIPSTask('CVEL')
    cvel.indata = indata
    cvel.source[1:] = cvelsource
    cvel.outna = cvel.inna
    cvel.outcl = cvel.incl
    cvel.gainuse = 7
    cvel.freqid = 1
    cvel.outseq = cvel.inseq+1
    cvel.outdisk = cvel.indisk
    cvel.doband  = doband
    cvel.bpver   = bpver
    cvel.aparm[10] = 1 # "trust cvel"

    cveldata = AIPSUVData(cvel.inna, cvel.incl, int(cvel.indisk), 2)
    if cveldata.exists():
        cveldata.clrstat()
        cveldata.zap()
    cvel()

    #uvsort(cveldata)
    indxr = AIPSTask('indxr')
    indxr.indata = cveldata
    indxr()

##############################################################################
# Manual phase cal
##############################################################################

def man_pcal(indata, refant, mp_source, mp_timera, debug, logfile, dpfour):

    if mp_source == ['']:
        mp_source = []
        for source in indata.sources:
            if source[0]=='F':
                mp_source.append(source)
    fringe            = AIPSTask('FRING')
    fringe.indata     = indata
    fringe.refant     = refant
    fringe.docal      = 1
    fringe.solint     = 6
    fringe.bchan      = 0
    fringe.echan      = 0
    fringe.aparm[1:]  =[2,0]
    fringe.dparm[2]   = 250
    fringe.dparm[3]   = 50
    fringe.dparm[4]   = dpfour
    fringe.dparm[8]   = 1
    fringe.snver      = 0
    #fringe.calso[1:]  = mp_source
    #fringe.inputs()
    #if mp_timera==0:
    #    sys.exit('No manual phase cal time range given!')
    #else:
    #    fringe.timer[1:] = mp_timera
    #    fringe()

    #qualfile = indata.name+'.'+indata.klass+'-qual.dat'
    #(source,timerange)=get_best_scan(indata,logfile, qualfile, 1)

    # Delete SN table from test fringe.
    #check_sncl(indata, 2, 6, logfile)
    fringe.calsour[1]  = mp_source
    fringe.timerang[1:] = mp_timerange
    fringe()

    sn=indata.table('AIPS SN', 0)
    mprint('###########################################',logfile)
    mprint('Found solutions for '+str(len(sn))+' of '
           +str(len(indata.antennas))+' antennas.',logfile)
    mprint('###########################################',logfile)
    return source, timerange

##############################################################################

def mprint(intext,logfile):
    print intext
    f=open(logfile, 'a')
    f.writelines(intext+'\n')
    f.close()

##############################################################################
# Table saving and restoring
##############################################################################

def runtacop(indata, outdata, inext, inver, outver, ncount):
    tacop         = AIPSTask('TACOP')
    tacop.indata  = indata
    tacop.outdata = outdata
    tacop.inext   = inext
    tacop.inver   = inver
    tacop.outver  = outver
    tacop.ncount  = ncount
    tacop()

def restore_su(indata, logfile):
    tasav_data=AIPSUVData(indata.name,'TASAV'+str(i),int(indata.disk),1)
    if tasav_data.exists():
        mprint('##########################################', logfile)
        mprint('TASAV file exists, restoring old SU table.', logfile)
        mprint('##########################################', logfile)
        while indata.table_highver('AIPS SU')>0:
            indata.zap_table('AIPS SU', 1)
        runtacop(tasav_data, indata, 'SU', 1, 1, 0)
    else:
        mprint('###############################################',logfile)
        mprint('No TASAV file. Restoring SU table not possible.',logfile)
        mprint('###############################################',logfile)

def restore_fg(indata, logfile):
    tasav_data=AIPSUVData(indata.name,'TASAV'+str(i),int(indata.disk),1)
    if tasav_data.exists():
        mprint('##########################################', logfile)
        mprint('TASAV file exists, restoring old FG table.', logfile)
        mprint('##########################################', logfile)
        while indata.table_highver('AIPS FG')>0:
            indata.zap_table('AIPS FG', 0)
        runtacop(tasav_data, indata, 'FG', 1, 1, 0)
    else:
        mprint('###############################################',logfile)
        mprint('No TASAV file. Restoring SU table not possible.',logfile)
        mprint('###############################################',logfile)

##############################################################################

def get_time():
    t=range(6)
    t[0]=time.localtime()[0]
    t[0]=str(t[0])
    for i in range(1,6):
        t[i]=time.localtime()[i]
        if t[i]<10:
            t[i]='0'+str(t[i])
        else:
            t[i]=str(t[i])
    a=t[3]+':'+t[4]+':'+t[5]+' on '+t[0]+'/'+t[1]+'/'+t[2]
    return a

##############################################################################
# Position shifting
##############################################################################

def shift_pos(indata, source, ra, dec, inver, outver):
    if source == '':
        source=findcal(indata, '')
    clcor             = AIPSTask('CLCOR')
    clcor.indata      = indata
    clcor.source[1]   = source
    clcor.opcode      = 'ANTC'
    clcor.clcorprm[5] = ra
    clcor.clcorprm[6] = dec
    clcor.gainver     = inver
    clcor.gainuse     = outver
    if ra!=0 or dec!=0:
        clcor()

##############################################################################

def runclcal(indata,snver,gainver,gainuse,interpol,dobtween,refant):
    clcal          = AIPSTask('CLCAL')
    clcal.indata   = indata
    clcal.refant   = refant
    clcal.snver    = snver
    clcal.inver    = 0
    clcal.gainver  = gainver
    clcal.gainuse  = gainuse
    clcal.interpol = interpol
    clcal.dobtween = dobtween
    clcal()

