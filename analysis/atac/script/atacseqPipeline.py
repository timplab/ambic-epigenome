#!/usr/bin/python3.5
'''
this is the pipeline for atac-seq
for now just making it applicable for 1x75 HiSeq (specifically for ambic)
e.g. ./atacseqPipeline.py -k /path/to/key.csv -o /path/to/working/dir
'''
import sys
import os
import subprocess

def parseArgs():
    import argparse
    ##first get dir of source code
    srcpath=sys.argv[0]
    srcdir=os.path.dirname(os.path.abspath(srcpath))
    ##parse arguments
    parser=argparse.ArgumentParser(description="atac-seq pipeline")
    parser.add_argument('-k','--key',type=os.path.abspath,required=True,help='path to the key csv; csv fields should be: 1)filepath,2)sampleID,3)replicate(group),4)(optional)read')
    parser.add_argument('-o','--output',type=os.path.abspath,help='directory for all outputs',required=True)
    parser.add_argument('-x','--btidx',type=str,help='bowtie2 index',required=True)
    parser.add_argument('-g','--genome',type=str,help='genome size for mac',required=True)
    parser.add_argument('-m','--module',type=str,default="all",help='run only certain module - for debugging, options:all,trim,align,peak')
    args=parser.parse_args()
    args.srcdir=srcdir
    # change the working directory
    os.chdir(args.output)
    print("Working directory is {}".format(os.getcwd()))
    # load csv and get info
    return args

class atacExperiment:
    ''' this parses the key file to determine the experiment parameters 
    inputs are 1) path of the key, 2) workinng directory, 3) source code directory'''
    def __init__(self,arg):
        self.btidx=arg.btidx
        self.srcdir=arg.srcdir
        self.wdir=arg.output
        self.keypath=arg.key
        self.genome=arg.genome
        self.content=self.readKey()
        self.spath,self.group,self.rep,self.base=self.parseKey()
        self.groups=list(set(self.group)) # distinct groups
        self.gnum=len(self.group) # number of distinct groups
        self.reps=[self.group.count(s) for s in self.groups] # number of replicates per group
        self.groupbase=[" ".join(x) for x in [[self.base[i] for i in range(0,len(self.base)) if x in self.group[i]]for x in self.groups]]
    def readKey(self):
        with open(self.keypath,'r') as fh:
#            fh.readline()
            data = [s.split(",") for s in [l.strip() for l in fh.readlines()]]
        return data
    def parseKey(self):
        fpath=[f[0] for f in self.content]
        group=[f[1] for f in self.content]
        rep=[f[2] for f in self.content]
        base=[os.path.basename(f[0]).split(".")[0] for f in self.content]
        return fpath,group,rep,base
    def isPE(self):
        ''' determine whether this run is paired-end sequenced based on input files'''
        try : 
            self.content[1][3]
            return True
        except IndexError:
            return False

def runparallel(runlist):
    p=[subprocess.Popen(fields) for fields in runlist]
    [x.wait() for x in p]

def runprocess(runlist):
    for fields in runlist:
        p = subprocess.Popen(fields)
        p.wait()

def trimAtac(run):
    outdir=os.path.join(run.wdir,"trim")
    fqdir=os.path.dirname(run.spath[1])
    if not os.path.exists(outdir):os.mkdir(outdir)
    trimpath=os.path.join(run.srcdir,"atacTrim.sh")
    trimarg=[trimpath,"-o",outdir,"-d",fqdir,"-s"]
    trimruns=[trimarg+[x] for x in run.base]
    runparallel(trimruns)
    print("done with trimming!")

def alignAtac(run):
    trimdir=os.path.join(run.wdir,"trim")
    outdir=os.path.join(run.wdir,"bam")
    if not os.path.exists(outdir):os.mkdir(outdir)
    alignpath=os.path.join(run.srcdir,"atacAlign.sh")
    alignarg=[alignpath,"-o",outdir,"-i",trimdir,"-x",run.btidx,"-t","10"]
    alignruns=[alignarg+["-s",x] for x in run.base]
    runprocess(alignruns)
    print("done with alignment")

def bamProcess(run):
    bamdir=os.path.join(run.wdir,"bam")
    processpath=os.path.join(run.srcdir,"atacbamProcess.sh")
    processarg=[processpath,"-o",bamdir]
    processruns=[processarg+["-s",x] for x in run.base]
    runprocess(processruns)
    # now generation beds
    bedpath=os.path.join(run.srcdir,"atacbedGeneration.sh")
    beddir=os.path.join(run.wdir,"bed")
    if not os.path.exists(beddir):os.mkdir(beddir)
    bedarg=[bedpath,"-i",bamdir,"-o",beddir]
    bedruns=[bedarg+["-g",x] for x in run.groups]
    runprocess(bedruns)
    print("done with processing bam files")

def peakAtac(run):
    beddir=os.path.join(run.wdir,"bed")
    outdir=os.path.join(run.wdir,"peak")
    if not os.path.exists(outdir):os.mkdir(outdir)
    peakpath=os.path.join(run.srcdir,"atacPeak.sh")
    peakarg=[peakpath,"-i",beddir,"-o",outdir,"-g",run.genome]
    peakruns=[peakarg+["-s",x] for x in run.base]
    runprocess(peakruns)
    print("done with peak calling")

def atacIDR(run):
    peakdir=os.path.join(run.wdir,"peak")
    outdir=os.path.join(run.wdir,"idr")
    if not os.path.exists(outdir):os.mkdir(outdir)
    idrpath=os.path.join(run.srcdir,"atacIDR.sh")
    idrarg=[idrpath,"-i",peakdir,"-o",outdir,"-b"]
    idrruns=[idrarg+[x] for x in run.groupbase]
    runprocess(idrruns)
    print("done with IDR")
    
    
def processAtac(run,module):
    if module=="all":
        module=['trim','align','bamprocess','peak','idr']
    if 'trim' in module:
        print("trimming")
        trimAtac(run)
    if 'align' in module:
        print("alignment")
        alignAtac(run)
    if 'bamprocess' in module:
        print("processing bam files")
        bamProcess(run)
    if 'peak' in module:
        print("peak calling")
        peakAtac(run)
    if 'idr' in module:
        print("performing IDR")
        atacIDR(run)


def main():
    args=parseArgs()
#    print(args)
    atacrun=atacExperiment(args)
    processAtac(atacrun,args.module)
    return 

if __name__ == "__main__":
    main()
