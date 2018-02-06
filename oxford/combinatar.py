#!/usr/bin/python

# using python 2.7.12
# usage e.g.: ./combineRuns.py [options] -o /path/to/output /path/to/runs
# options:
#   -n number of files per dir
#   -s name of the run to combine
#   -o output directory
#   -a archive

def parseArgs():
    import argparse
    import os
    parser=argparse.ArgumentParser(description="Combine several nanopore runs into one dir")
    parser.add_argument('root',metavar='dir',type=os.path.abspath,help='common directory for fast5 files')
    parser.add_argument('-n',metavar='N',type=int,default=4000,help='number of fast5s to group into one dir')
    parser.add_argument('-s','--sample',metavar='foo',type=str,required=False,help='name of the run to combine')
    parser.add_argument('-a','--archive',action='store_true',help='flag for archiving')
    parser.add_argument('-o','--outdir',metavar='dir',type=os.path.abspath,required=True,help='output directory')
    args=parser.parse_args()
    if args.sample==None:
        ## if sample name is not given, using the dir above root and longest suffix of the dirs as the prefix
        dirs=os.listdir(args.root)
        dirsrev=[s[::-1] for s in dirs]
        lcsrev=os.path.commonprefix(dirsrev)
        lcs=lcsrev[::-1].strip("_")
        args.sample=lcs
        print("sample name not given, using {} as sample name".format(lcs))
    if not os.path.exists(os.path.join(args.outdir,args.sample)):
        os.mkdir(os.path.join(args.outdir,args.sample))
    return args

def indUpdate(wdir,current,last):
    import os
    current+=1
    if current==last:
        i=str(int(os.path.basename(wdir))+1)
        root=os.path.dirname(wdir)
        wdir=os.path.join(root,i)
        current=0
        os.mkdir(wdir)
    return wdir,current

def combineRuns(args):
    import os
    n=args.n
    outdir=os.path.join(args.outdir,args.sample) # confusion, but this new outdir contains sample dir
    #i = output folder index, j = output file index, x = total file count
    if os.path.exists(outdir):
        dirs=os.listdir(outdir)
        folder_ind=len(dirs)
    else:
        folder_ind= 0
    ind=str(folder_ind)
    wdir=os.path.join(outdir,ind)
    os.mkdir(wdir)
    file_ind = -1
    totfiles=0
    for root,dirs,files in os.walk(args.root):
        for f in files:
            if args.sample in f:
                wdir,file_ind=indUpdate(wdir,file_ind,n)
                path=os.path.join(root,f)
                newpath=os.path.join(wdir,f)
                os.rename(path,newpath)
                totfiles+=1
    print("total number of files: {}".format(totfiles))
    if args.archive:
        import tarfile
        print("archiving...")
        with tarfile.open(outdir+".raw.tgz","w:gz") as fh:
            for i in next(os.walk(outdir))[1]:
                print(i)
                fh.add(os.path.join(outdir,i),arcname=i)

def main():
    args=parseArgs()
    combineRuns(args)

if __name__ == "__main__":
    main()
