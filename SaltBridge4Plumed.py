#/usr/bin/env python

import sys
import os
# this script generates a plumed input file for salt bridge determination
class PdbIndex :
    """
    Input the reisude number sequence, then out put the required index atoms
    """

    def __init__(self, inpdb, chain, atomtype, residueNdx) :
        self.chain = chain
        self.inpdb = inpdb
        self.atomtype = atomtype
        self.residueNdx = residueNdx

    def res_index(self, atomList, atomtype="", dihedraltype="None"):

        #atomInfor = {}
        indexlist = []

        if atomtype == "dihedral" :
            # format of indexlist : [ [1, 5, 7, 9 ], ]
            # residueNdx = [ 5, 6, 7, 8 ]
            #print "dihedral here"
            indexlist = []
            print residueNdx
            for resndx in residueNdx :
                #print "PHI here"
                #define index for phi angle C N CA C
                phitype = ['C', 'N', 'CA', 'C']
                phipair = [-1, -1, -1, -1]
                phindxadd = [-1, 0, 0, 0]
                for i in range(4) :
                    with open(self.inpdb) as lines :
                        for s in lines :
                            if len(s.split()) > 1 and s.split()[0] == "ATOM" and s[21] == self.chain :
                                if int(s[22:26].strip()) == resndx + phindxadd[i] and s[12:16].strip() == phitype[i] :
                                    phipair[i] = int(s.split()[1])
                                    #print int(s.split()[1])
                                    #break

                psitype = ['N', 'CA', 'C', 'N']
                psindxadd = [ 0, 0, 0, 1]
                psipair = [-1, -1, -1, -1]
                for i in range(4):
                    with open(self.inpdb) as lines:
                        for s in lines:
                            if len(s.split()) > 1 and s.split()[0] == "ATOM" and s[21] == self.chain:
                                if int(s[22:26].strip()) == resndx + psindxadd[i] and s[12:16].strip() == psitype[i]:
                                    psipair[i] = int(s.split()[1])
                                    #break

                if "PHI" in dihedraltype :
                    indexlist.append(psipair)
                if "PSI" in dihedraltype :
                    indexlist.append(phipair)

        else :

            with open(self.inpdb) as lines :
                for s in lines :
                    if len(s.split()) > 1 and s.split()[0] == "ATOM" and s[21] == self.chain :
                        if int(s[22:26].strip()) in self.residueNdx :
                            if atomtype == "non-hydrogen" :
                                if s[13] != "H" and s.split()[2][0] != "H" and s.split()[-1] != "H" :
                                    indexlist.append(s.split()[1])
                            elif atomtype == "all-atom" :
                                ## all atoms
                                indexlist.append(s.split()[1])
                            elif atomtype == "side-chain" :
                                if s[12:16].strip() not in ['CA', 'N', 'C', 'O'] :
                                    indexlist.append(s.split()[1])
                            else:
                                if s[12:16].strip() in atomList :
                                    indexlist.append(s.split()[1])
                        else :
                            pass
                    else:
                        pass

        return indexlist

    def atomList(self):
        atomList = []
        atomtype = ""
        if "mainchain" in self.atomtype or "Main" in self.atomtype or "main" in self.atomtype :
            atomList = ['CA', 'N', 'C', 'O']

        elif "CA" in self.atomtype or "ca" in self.atomtype or "Ca" in self.atomtype or "alpha" in self.atomtype :
            atomList = ['CA']

        elif "backbone" in self.atomtype or "Back" in self.atomtype or "bone" in self.atomtype :
            atomList = ['CA', 'N']

        elif "all" in self.atomtype :
            atomtype = "all-atom"

        elif "no hy" in self.atomtype or "non-hy" in self.atomtype :
            atomtype = "non-hydrogen"

        elif "side" in self.atomtype or "Sidechain" in self.atomtype or "sidechain" in self.atomtype :
            atomtype = "side-chain"

        elif "PSI" in self.atomtype or "PHI" in self.atomtype or "phi" in self.atomtype or 'psi' in self.atomtype :
            atomtype = "dihedral"

        return  atomList, atomtype


class SaltBridgeResidues :
    # define the atom index of the salt bridges
    def __init__(self, input):
        self.reference = input

    def residueNdx(self, resNameList=['ALA'], resNdxRange=[0,0], chains=['A']):
        resIndex =[]

        with open(self.reference) as lines :
            for s in lines :
                if ("ATOM" in s) and s[21] in chains \
                        and int(s[22:26].strip()) in resNdxRange \
                        and s[17:20].strip() in resNameList :
                    resIndex.append(int(s[22:26].strip()))
        return resIndex

    def saltbridgePlumedInput(self, atomSet1, atomSet2,
                           distCoutoff=0.32,
                           poutFileName="saltbridge.dat",
                           sboutFileName="num_sb_vs_time.dat",
                           ):
        with open(poutFileName, 'wb') as pf :
            pf.write("# A plumed input for counting number of salt-bridges \n")

            # define groups
            pf.write("groupA: GROUP ATOMS=%s" % atomSet1[0])
            for atom in atomSet1[1:] :
                pf.write(",%s"%atom)
            pf.write(" \n")

            # define groups
            pf.write("groupB: GROUP ATOMS=%s"% atomSet2[0])
            for atom in atomSet2[1:] :
                pf.write(",%s"%atom)
            pf.write(" \n")

            # define coordination number here
            pf.write("sb: COORDINATION GROUPA=groupA GROUPB=groupB "
                     "R_0=%4.2f NLIST NL_CUTOFF=%4.2f NL_STRIDE=5 \n" %
                     (distCoutoff, 2.0*distCoutoff)
                     )

            # output data
            pf.write("PRINT ARG=sb STRIDE=1 FILE=%s \n"%sboutFileName)

if __name__ == "__main__" :
    # cd $PWD
    os.chdir(os.getcwd())

    args = sys.argv
    if len(args) <= 1 :
        d = '''
        Generating Input For Plumed to Calculate Number of Salt Bridges vs Time.
        Usage:

        python SaltBridge4Plumed.py reference.pdb plumed.dat
        '''
        print(d)
        sys.exit(1)

    pdbfile = args[1]

    # group A
    ares = SaltBridgeResidues(pdbfile)
    resndx = ares.residueNdx(['ARG','LYS',], range(4, 1400), ['B'])
    andx = PdbIndex(pdbfile, "B", "", resndx)
    aatom = andx.res_index(['NH1'], )

    # group B
    bres = SaltBridgeResidues(pdbfile)
    batom = []
    for c in ['A','C','D'] :
        resndx = bres.residueNdx(['DA','DT','DC','DG','A','U','C','G'], range(1, 200), [c])
        bndx = PdbIndex(pdbfile, c, "", resndx)
        batom += bndx.res_index(['O1P'])

    cplxsb = SaltBridgeResidues(pdbfile)
    cplxsb.saltbridgePlumedInput(aatom, batom, 0.32, 'plumed_sb.dat')

    print("Run plumed:\nplumed driver --mf_xtc your.xtc --plumed plumed.dat")