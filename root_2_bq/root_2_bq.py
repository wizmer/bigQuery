import os
import sys
import glob
import gzip
import pandas as pd
import ROOT
import root_numpy

workdir = "/afs/cern.ch/user/k/kostams/eos/ams/user/k/kostams/DeuteronNtuples/production/"
filenames = glob.glob(workdir + "*.ntuple.root")

def steps(N,n):
    current = 0
    while current < N:
        yield current, current+n
        current += n

for n,(i,j) in enumerate(steps(len(filenames),10)):
    
    outfname = "fulldata." + str(n) + ".csv.gz" 
    print "Current target is " + outfname
    if os.path.exists(workdir + outfname):
        print " already exists skipping."
        continue
    os.system("touch " + workdir + outfname)
    
    data = []
    for filename in filenames[i:j]:
        print "Reading " + filename.split('/')[-1]
        sys.stdout.flush()
        f = ROOT.TFile(filename)
        if f.IsZombie(): continue
        tree = f.Get("selections")
        data.append(pd.DataFrame(root_numpy.tree2rec(tree)))
    data = pd.concat(data)

    print "Writing " + str(len(data)) + " rows into " + outfname + " ... ",
    sys.stdout.flush()
    with gzip.open(outfname, 'w') as csvOut: data.to_csv(csvOut)  
    print " done."
    
    print "Copying ...",
    sys.stdout.flush()
    os.system("cp " + outfname + " " + workdir)
    print " done."
    
    print "Cleanup ...",
    sys.stdout.flush()
    os.system("rm " + outfname)
    print " done."
