#/usr/bin/env python

import sys
import os
import argparse
from argparse import RawTextHelpFormatter

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

def arguments() :
    d = '''
    Prepare plumed input file for salt-bridge determination.
    '''
    parser = argparse.ArgumentParser(description=d)

    parser.add_argument('-ref', '--reference', type=str, default='reference.pdb',
                        help='Reference PDB file to get index of atoms. \n'
                             'Default is reference.pdb ')

    parser.add_argument('-mt','--moleculetype', type=str, nargs='+', default=['Protein', 'Protein'],
                        help="Select molecule types for your two molecules. \n"
                             "Default is: Protein Protein. ")
    parser.add_argument('-rn','--residuename', type=str, nargs='+', default=['ARG','LYS'],
                        help="By default, if a protein involves salt-bridges, we choose "
                             "LYS and ARG for distance calculating. \nYou could specify your own.")
    parser.add_argument('-o', '--output', default='plumed.dat',type=str,
                        help="Output file name for plumed. \n Default is Plumed.dat.")
    parser.add_argument('-res', '--residuerange',default=[], nargs='+',
                        help="Residue sequence index range for residues involving salt bridges. ")
    parser.add_argument('-cA', '--chainsA', type=str, default=['A'], nargs='+',
                        help="For molecule A, select the chains of molecule A.")
    parser.add_argument('-cB', '--chainsB', type=str, default=['B'], nargs='+',
                        help="For molecule B, select the chains of molecule B.")
    parser.add_argument('-cutoff','--distance_cutoff', default=0.32, type=float,
                        help="Distance cutoff (nanometer) for an effective salt-bridge like contact. \n"
                             "Default is 0.32 nm. ")

    args = parser.parse_known_args()

    return(args)

if __name__ == "__main__" :
    # cd $PWD
    os.chdir(os.getcwd())

    args = sys.argv
    if len(args) <= 1 :
        d = '''
        Generating Input For Plumed to Calculate Number of Salt Bridges vs Time.
        Usage:

        python SaltBridge4Plumed.py -h
        '''
        print(d)
        sys.exit(1)

    args = arguments()

    pdbfile = args.ref

    # group A B
    atomsAB = []
    try:

        for i in range(2) :
            aatom = []

            chains = [args.cA, args.cB]

            ares = SaltBridgeResidues(pdbfile)
            resndx = ares.residueNdx(args.rn, range(args.res[0], args.res[1]), args.cA)

            for chain in chains[i] :
                    andx = PdbIndex(pdbfile, chain, "", resndx)
                    if args.mt[i] in ["Protein", "protein", "pro"]:
                        if i == 0 :
                            aatom += andx.res_index(['NH1', 'NH2', 'NE'], )
                        else :
                            aatom += andx.res_index(['O',], )
                    elif args.mt[i] in ['RNA', 'DNA', 'rna', 'dna', 'nucleic'] :
                        aatom += andx.res_index(['OP1','OP2'], )
            atomsAB.append(aatom)
    except :
        print("Name Error!")

    cplxsb = SaltBridgeResidues(pdbfile)
    cplxsb.saltbridgePlumedInput(atomsAB[0], atomsAB[1], args.cutoff, args.o)

    print("Run plumed:\nplumed driver --mf_xtc your.xtc --plumed %s\n"%args.o)