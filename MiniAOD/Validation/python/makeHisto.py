import ROOT
import sys, array, math

from optparse import OptionParser

usage = "Usage: %prog [options] sample1.root sample2.root"
parser = OptionParser(usage=usage)
parser.add_option("-p","--plot",        default="plot.def", help="File name with plot definition. \nDefault: %default", action="store_true")
parser.add_option("-o","--outFile",     default="output.root", help="Output file to be produced. \nDefault: %default")
(options,arg) = parser.parse_args()

args = sys.argv[1:]

if len(args) != 2:
    parser.error("incorrect number of arguments")

branchNames = []
branchList = []
histos = []

def DeltaR(e1, e2, p1, p2):
    dp = abs(p1-p2)
    if (dp > math.pi):
        dp = dp - (2*math.pi)
    return math.sqrt((e1-e2)*(e1-e2) + dp*dp)

def listOfBranches(t, n):
    global branchNames, branchList, histos, options

    desc = open(options.plot)
    lines = desc.readlines()
    desc.close()
    
    bntemp = []
    bltemp = {}
    htemp = {}

    for b in t.GetListOfLeaves():
        bname = b.GetName()
        for l in lines:
            fields = l.split()
            if (bname == fields[0]):
                bType = t.GetLeaf(bname).GetTypeName()
                if (bType == "Float_t"):
                    bltemp[bname] = array.array('f', [0,0,0,0,0,0,0,0,0,0])
                elif (bType == "Int_t"):
                    bltemp[bname] = array.array('i', [0,0,0,0,0,0,0,0,0,0])
                elif (bType == "Double_t"):
                    bltemp[bname] = array.array('d', [0,0,0,0,0,0,0,0,0,0])
                else:
                    continue
                bntemp.append(bname)
                t.SetBranchAddress(bntemp[-1], bltemp[bntemp[-1]])
                htemp[bntemp[-1]] = ROOT.TH1F(bntemp[-1]+str(n), "", int(fields[1]), float(fields[2]), float(fields[3]))

    branchNames.append(bntemp)
    branchList.append(bltemp)
    histos.append(htemp)

fnames = []
fnames.append(args[0])
fnames.append(args[1])
files = []
trees = []
m = []

for fname in fnames:
    mymap = {}
    files.append(ROOT.TFile(fname))
    trees.append(files[-1].Get("validation"))
    entries = trees[-1].GetEntries()

    for z in xrange(entries):
        trees[-1].GetEntry(z)
        mymap[(trees[-1].run, trees[-1].event)] = z
    m.append(mymap)

for i in xrange(2):
    listOfBranches(trees[i], i)

for k in m[0].keys():
    trees[0].GetEntry(m[0][k])
    trees[1].GetEntry(m[1][k])
    for n1 in xrange(trees[0].n):
        eta1 = trees[0].eta[n1]
        phi1 = trees[0].phi[n1]
        for n2 in xrange(trees[1].n):
            eta2 = trees[1].eta[n2]
            phi2 = trees[1].phi[n2]
            
            dR = DeltaR(eta1, eta2, phi1, phi2)
            if (dR < 0.1):
                for b in branchNames[0]:
                    histos[0][b].Fill(branchList[0][b][n1])
                    histos[1][b].Fill(branchList[1][b][n2])
                break

out = ROOT.TFile(options.outFile, "recreate")
for i in xrange(2):
    for k in histos[i].keys():
        histos[i][k].Write()
out.Close()
        
    


