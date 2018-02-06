#! /usr/bin/env python3

class h5read:
    
    def __init__(self, filename):
        self.fastq=[]
        self.stime=-1
        self.etime=-1
        self.duration=-1
        self.rname=[]
        self.bseq=[]
        self.qual=[]
        self.avgqual=-1
        self.rlength=-1
        self.chnum=-1
        self.mux=-1
        self.sampling=-1
        self.fname=filename
        self.bc=-1

        ##Setup for 1D only

        self.fqkey='/Analyses/Basecall_1D_000/BaseCalled_template/Fastq'
        
        self.hdf = h5py.File(filename, 'r')
        
    def check(self):
        ##Check to see if things are in hdf file
        if self.fqkey in self.hdf:
            self.fqcheck=True
        else:
            self.fqcheck=False

    def fqextract(self):
        import statistics

                
        self.fastq = self.hdf[self.fqkey][()]
        ##Python3 wants explicit for ascii, assumed unicode
        fqlines=self.fastq.decode('ascii').split('\n')

        self.rname=fqlines[0]
        self.bseq=fqlines[1]
        self.qual=fqlines[3]

        self.rlength=len(self.qual)

        qualval = [ord(x) - 33 for x in self.qual.rstrip('\n')]
        self.avgqual = statistics.mean(qualval)

    def read_rawblock_get(self):
        ##Shamelessly stolen/adapted from poretools source
        raw_reads = self.hdf.get('Raw/Reads')

        reads = list(raw_reads.keys())
        if len(reads)==0:
            print("No reads in raw")
            return None
        elif len(reads)>1:
            print("More than one read?!")
            return None
        else:
            readpath= 'Raw/Reads/%s' % ( reads[0] )
            return self.hdf.get(readpath)

    def textract(self):
        ##Note - only works on relatively recent data

        globalkey=self.hdf['/UniqueGlobalKey']

        readblock=self.read_rawblock_get()
        
        ##get times this way
        expstarttime=globalkey['tracking_id/'].attrs.get('exp_start_time').decode('ascii')

        if expstarttime.endswith('Z'):
            ##Parse this datetime
            expparsed=dateutil.parser.parse(expstarttime)
            self.expstime=int(time.mktime(expparsed.timetuple()))
            
            self.sampling=globalkey['channel_id'].attrs['sampling_rate']
            
            self.duration=readblock.attrs['duration'] / self.sampling
            
            self.stime=readblock.attrs['start_time'] / self.sampling + self.expstime
            
            self.chnum=globalkey['channel_id'].attrs['channel_number'].decode('ascii')
            self.mux=readblock.attrs['start_mux']
            
            self.etime=self.stime+self.duration
        

    def bcextract(self):

        bcsummary=self.hdf['/Analyses/Barcoding_000/Summary']
        bcstring=bcsummary['barcoding'].attrs['barcode_arrangement']

        #print(bcstring)
        #print(type(bcstring))
        if (bcstring=="unclassified"):
            self.bc=0
        else:
            self.bc=int(bcstring[7:])
            
    def fqout(self, fout):
        rnameout=self.rname.split(' ')[0]+":1D_000:template "+os.path.basename(self.fname)+" "+self.fname
        #print(rnameout)
        fout.write(self.fastq)


    def timeout(self, tout):
        delim=','
        tout.write(delim.join([self.rname,str(self.rlength),str(self.avgqual), str(self.stime), str(self.etime), str(self.duration), str(self.chnum), str(self.mux), str(self.bc)])+'\n')



class seqout:

    def __init__(self, prefix, barcode):
        self.opened=False
        ## -1 barcode means no barcode
        ## 0 barcode is unclassified
        self.barcode=barcode
        if (barcode < 0):
            self.name=prefix+'.fq.gz'
        else:
            self.name=prefix+'.'+str(barcode)+'.fq.gz'
        
        self.fout=[]

    def fqopen(self):
        if not (self.opened):
            self.fout=gzip.open(self.name, 'wb')
            self.opened=True

    def fqclose(self):
        if (self.opened):
            self.fout.close
            self.opened=False


##Main 

##HDF format
import h5py
##Os IO stuff
import os
import shutil
##Tarball 
import tarfile
##File recognize
import glob
##shell util
import shutil
##Random number
import random
##Gzip lib
import gzip
##Stats
import statistics
##Get args
import argparse
##RE
import re
##Tempfile directory generate
import tempfile
##For timestamp manip
import dateutil.parser
import time

parser = argparse.ArgumentParser( description='Extract fastq from fast5 as well as times')
parser.add_argument('--input', '-i', type=str, required=True, help='input location of folder fast5 files')
parser.add_argument('--barcode', '-b',  action='store_true', help='Sample is barcoded')
parser.add_argument('--time', '-t', action='store_true', help='Put out times')
parser.add_argument('--output', '-o', type=str, required=True, help='output filename prefix')


args=parser.parse_args()

filelist=glob.glob(args.input+'/*.fast5', recursive=True)


##prefix, barcode
if (args.barcode):
    fqoutfile=[seqout(args.output, i) for i in range(97)]
else:
    fqoutfile=seqout(args.output, -1)
    fqoutfile.fqopen()

if args.time:
    tout=gzip.open(args.output+'.tout.csv.gz', 'wt')
else:
    tout=''


for filename in filelist:
    try: 
        ##Load file
        read=h5read(filename)
        read.check()

        if (read.fqcheck):
            read.fqextract()
            ##Get barcode
            ##Need to write somehow specific to if barcode
            if (args.barcode):
                read.bcextract()
                fqoutfile[read.bc].fqopen()
                read.fqout(fqoutfile[read.bc].fout)
            else:
                read.fqout(fqoutfile.fout)
            if (args.time):
                read.textract()
            if (args.time):
                read.timeout(tout)            
        read.hdf.close()        
    except Exception as e:
        print(e)
        print(filename+" failed to open!!")
        
if (args.barcode):
    for i in range(97):
        fqoutfile[i].fqclose()
else:
    fqoutfile.fqclose()
    

if (args.time):
    tout.close()
