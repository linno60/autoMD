from autoMD import *

ligname = "AS0"

spdb = SummaryPDB('2ovh.pdb','amino-acid.lib')
prot, misndx, noprondx, noproname, fullndx = spdb.summary('A', verbose=True)

print misndx
rec_ndx = [int(misndx['A'][0]), int(misndx['A'][-1])]
print noproname['A']
for item in noproname['A'] :
    if ligname in item :
        ligNdx = int(item.split("_")[-1])

print ligNdx
print rec_ndx

fs, misslist, fullndx = spdb.missingRes('A', '', verbose=True, pdbCode='2ovh')
print misslist

fpdb = FixPDB()
file = fpdb.addLoopLOOPY('2ovh', rec_ndx, misslist.values()[0], 'A',
                         loopyexe='/home/liangzhen/software/loopy/test/loopy',ligandchain='A',
                         ligandName='AS0',ligandNdx=ligNdx, verbose=True)
print "Best Model File name is"
print file