import ROOT
import sys

from optparse import OptionParser

usage = "Usage: %prog [options] plotfile.root"
parser = OptionParser(usage=usage)
parser.add_option("-p","--plot",        default="plot.def", help="File name with plot definition. \nDefault: %default", action="store_true")
parser.add_option("-o","--outFile",     default="plots.root", help="Output file to be produced. \nDefault: %default")
(options,arg) = parser.parse_args()

args = sys.argv[1:]

if len(args) != 1:
    parser.error("incorrect number of arguments")

ROOT.gROOT.SetBatch(True)
canvases = []
histos = [{}, {}]
plotnames = []

desc = open(options.plot)
lines = desc.readlines()
desc.close()
    
for l in lines:
    plotnames.append(l.split()[0])

f = ROOT.TFile(args[0])
for n in plotnames:
    histos[0][n] = f.Get(n+"0").Clone()
    histos[1][n] = f.Get(n+"1").Clone()

for k in histos[0]:
    canvases.append(ROOT.TCanvas("c" + k, "c"))
    histos[0][k].Draw()
    histos[1][k].Draw("PESAME")
    histos[0][k].SetLineColor(ROOT.kBlue)
    histos[1][k].SetLineColor(ROOT.kRed)
    histos[1][k].SetMarkerColor(ROOT.kRed)
    histos[1][k].SetMarkerStyle(21)
    histos[0][k].Draw("SAME")
    canvases[-1].SaveAs(options.outFile.replace(".root", "")+"_"+k+".png")

output = ROOT.TFile(options.outFile, "recreate")
for c in canvases:
    c.Write()
output.Close()
