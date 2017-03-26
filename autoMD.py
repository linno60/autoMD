#!/usr/bin/env python

import os
import sys
import urllib2
import time
import argparse
import subprocess as sp
from argparse import RawTextHelpFormatter
import math
from glob import glob
from collections import *
import linecache
import numpy as np

# import modeller for loop refinement
try:
    from modeller import *
    from modeller.automodel import *
    MODELLER_EXIST = True
except ImportError :
    print("Warning: Modeller is not well imported, some features may not work. ")
    MODELLER_EXIST = False

class PdbIndex :
    '''
    Input the reisude number sequence, then out put the required index atoms
    '''
    def __init__(self, ) :
        pass

    def res_index(self,inpdb, chain, atomtype, residueNdx, atomList, dihedraltype="None"):
        '''
        Obtain atom index from a reference pdb file
        provide information including: residue indexing, atomtype, atomname

        :param inpdb:
        :param chain:
        :param atomtype:
        :param residueNdx:
        :param atomList:
        :param dihedraltype:
        :param atomName:
        :return: a list, of atom index
        '''

        #atomInfor = {}
        indexlist = []
        if len(residueNdx) == 2 :
            residueNdx = range(residueNdx[0], residueNdx[1]+1)
        elif len(residueNdx) == 1 :
            residueNdx = range(residueNdx[0], residueNdx[0]+1)
        else :
            print("Error!! No residue index provided. ")
            residueNdx = []

        if atomtype == "dihedral" :
            indexlist = []
            for resndx in residueNdx :
                phitype = ['C', 'N', 'CA', 'C']
                phipair = [-1, -1, -1, -1]
                phindxadd = [-1, 0, 0, 0]
                for i in range(4) :
                    with open(inpdb) as lines :
                        for s in lines :
                            if len(s.split()) > 1 and s.split()[0] == "ATOM" and s[21] == chain :
                                if int(s[22:26].strip()) == resndx + phindxadd[i] and s[12:16].strip() == phitype[i] :
                                    phipair[i] = int(s.split()[1])

                psitype = ['N', 'CA', 'C', 'N']
                psindxadd = [ 0, 0, 0, 1]
                psipair = [-1, -1, -1, -1]
                for i in range(4):
                    with open(inpdb) as lines:
                        for s in lines:
                            if len(s.split()) > 1 and s.split()[0] == "ATOM" and s[21] == chain:
                                if int(s[22:26].strip()) == resndx + psindxadd[i] and s[12:16].strip() == psitype[i]:
                                    psipair[i] = int(s.split()[1])
                                    #break

                if "PHI" in dihedraltype :
                    indexlist.append(psipair)
                if "PSI" in dihedraltype :
                    indexlist.append(phipair)

        else :

            with open(inpdb) as lines :
                for s in lines :
                    if len(s.split()) > 1 and s.split()[0] == "ATOM" and s[21] == chain :
                        if int(s[22:26].strip()) in residueNdx :
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
#                                print "ATOM LIST"
                                if s[12:16].strip() in atomList :
                                    indexlist.append(s.split()[1])
                        else :
                            pass
                    else:
                        pass

        return(indexlist)

    def atomList(self, atomtype, atomname):
        '''
        provide information of atom type and atomname
        :param atomtype:
        :param atomname:
        :return:
        '''
        atomList = []
        if atomname :
            atomList = atomname
        else :
            if "mainchain" in atomtype or "Main" in atomtype or "main" in atomtype :
                atomList = ['CA', 'N', 'C', 'O']

            elif "CA" in atomtype or "ca" in atomtype or "Ca" in atomtype or "alpha" in atomtype :
                atomList = ['CA']

            elif "backbone" in atomtype or "Back" in atomtype or "bone" in atomtype :
                atomList = ['CA', 'N']

            elif "all" in atomtype :
                atomtype = "all-atom"

            elif "no hy" in atomtype or "non-hy" in atomtype :
                atomtype = "non-hydrogen"

            elif "side" in atomtype or "Sidechain" in atomtype or "sidechain" in atomtype :
                atomtype = "side-chain"

            elif "PSI" in atomtype or "PHI" in atomtype or "phi" in atomtype or 'psi' in atomtype :
                atomtype = "dihedral"

        return(atomList, atomtype)

    def atomList2File(self, atomNdxList, groupName, outputNdxFile, append=True):
        if append :
            tofile = open(outputNdxFile, 'a')
        else :
            tofile = open(outputNdxFile, 'wb')

        tofile.write('[ %s ] \n'%groupName)
        i = 0
        for atom in atomNdxList :
            i += 1
            tofile.write('%6d ' % int(atom))
            if i % 15 == 0:
                tofile.write('  \n')

        tofile.close()
        return(1)

    def arguements(self) :
        d = '''
        ################################################################
        # Generate GMX Index from a PDB file
        # Generate POSRES file for a PDB File or Gro File
        # Contact  LZHENG002@E.NTU.EDU.SG
        # Version  2.2
        # Update Mar 23, 2017
        ################################################################

        Usage examples

        Generate backbone index for residue 1 to 1000 within Chain B
        python autoMD.py index -res 1 100 -at backbone -chain B -out index.ndx

        Generate all atom index for residue 78 to 100
        python autoMD.py index -inp S_1.pdb -out index.ndx -at allatom -chain ' ' -res 78 100

        Generate dihedral (PHI only) quadroplex index
        python autoMD.py index -inp S_1.pdb -out index.ndx -at dihedral -chain ' ' -res 78 100 -dihedral PHI

        '''

        parser = argparse.ArgumentParser(description=d, formatter_class=RawTextHelpFormatter)
        #parser.add_argument('-h','--help', help="Show this help information. \n")
        parser.add_argument('-pdb', '--pdbfile',type=str, default='input.pdb',
                            help="Input PDB file for generating Index. ")
        parser.add_argument('-out', '--output',type=str, default='output.ndx',
                            help="Output Index file including atoms sequence number.\n"
                                 "Default name is output.ndx \n")
        parser.add_argument('-at', '--atomtype',type=str, default='allatom',
                            help="Selected atom type for generating index. \n"
                                 "Options including: allatom, mainchain, \n"
                                 "non-hydrogen, c-alpha, backbone, sidechain, dihedral\n"
                                 "Default choice is: allatom \n")
        parser.add_argument('-an','--atomname', type=str, default=[], nargs='+',
                            help="Select atom by atom names. A list of atom names \n"
                                 "could be supplied. If not given, all atom names are \n"
                                 "considered. Default is [].")
        parser.add_argument('-chain', '--chainId',type=str, default="A",
                            help="Protein chain identifier. Default chain ID is A. \n")
        parser.add_argument('-res', '--residueRange',type=int, nargs='+',
                            help="Residue sequence number for index generating. \n"
                                 "Example, -res 1 100, generateing atom index within \n"
                                 "residues 1 to 100. Default is None.")
        parser.add_argument('-posres', '--positionRestraint',default=False,
                            help="Generate a postion restraint file for selected index.\n"
                                 "Default name is posres.itp \n")
        parser.add_argument('-dihe', '--dihedralType',default=None, type=str,
                            help="Whether generate dihedral angles index (quadruplex).\n"
                                 "Phi and Psi are considered. Optional choices are: \n"
                                 "PHI, PSI, PHI_PSI, or NA. Default is NA. \n")
        parser.add_argument('-gn','--groupName', type=str, default=None,
                            help="The name of the group of atoms selected. \n"
                                 "Default is None.")
        parser.add_argument('-append', '--appendFile', default=True, type=bool,
                            help="Append the group of index to the output file. \n"
                                 "Options: True, False. \n"
                                 "Default is True.")

        args, unknown = parser.parse_known_args()

        # decide to print help message
        if len(sys.argv) < 3 :
            # no enough arguements, exit now
            parser.print_help()
            print("\nYou chose non of the arguement!\nDo nothing and exit now!\n")
            sys.exit(1)

        return(args)

    def genGMXIndex(self):
        '''
        run geneIndex
        initialize argument parser function.
        :return:
        '''
        args = self.arguements()

        if len(args.residueRange) == 2:
            residueNdx = range(args.residueRange[0], args.residueRange[-1] + 1)
        elif len(args.residueRange) == 1:
            residueNdx = range(args.residueRange[0], args.residueRange[0] + 1)
        elif len(args.residueRange) > 2:
            residueNdx = args.residueRange
        else:
            print "\nNumber of residue id is not correct. \nExit Now. "
            sys.exit(1)

        if args.dihedralType :

            atomlist, atomtype = self.atomList(args.dihedralType, args.atomname)
        else :
            atomlist, atomtype = self.atomList(args.atomtype, args.atomname)

        ## generate index
        atomndx = self.res_index(inpdb=args.pdbfile, chain=args.chainId,
                                 residueNdx= args.residueRange, atomList=atomlist,
                                 atomtype=atomtype, dihedraltype=args.dihedralType
                                 )
        if args.appendFile :
            append = 'a'
        else :
            append = 'wb'
        tofile = open(args.output, append)
        
        if args.groupName :
            tofile.write('[ %s ] \n' % (args.groupName.strip()))
        else :
            tofile.write('[ %s_%d_%d ] \n' % (args.chainId, args.residueRange[0], args.residueRange[-1]))

        if atomtype == "dihedral":
            for index in atomndx:
                for i in range(4):
                    tofile.write("%5d " % index[i])
                tofile.write("  \n")

        else:
            i = 0
            for atom in atomndx:
                i += 1
                tofile.write('%4s ' % atom)
                if i % 15 == 0:
                    tofile.write('  \n')
            tofile.write(" \n")
            tofile.close()

            if args.positionRestraint:
                with open('posres.itp', 'w') as tofile :
                    tofile.write("[ position_restraints ] \n; ai  funct  fcx    fcy    fcz  \n")
                    for atom in atomndx:
                        tofile.write("%12d  1  1000  1000  1000  \n" % int(atom))

        print("\nGenerating Index File Completed!")
        return(1)

class RewritePDB :
    def __init__(self, filename):
        self.pdb = filename

    def pdbRewrite(self, output, chain, atomStartNdx, resStartNdx):
        # change atom id, residue id and chain id
        resseq = resStartNdx
        atomseq = int(atomStartNdx)
        chainname = chain
        lines = open(self.pdb)
        newfile = open(output,'w')
        resseq_list = []

        for s in lines :
            if "ATOM" in s and len(s.split()) > 6 :
                atomseq += 1
                newline = s
                newline = self.atomSeqChanger(newline, atomseq)
                newline = self.chainIDChanger(newline, chainname)
                if len(resseq_list) == 0 :
                    newline = self.resSeqChanger(newline, resseq)
                    resseq_list.append(int(s[22:26].strip()))
                else :
                    if resseq_list[-1] == int(s[22:26].strip()) :
                        newline = self.resSeqChanger(newline, resseq)
                    else :
                        resseq += 1
                        newline = self.resSeqChanger(newline, resseq)
                    resseq_list.append(int(s[22:26].strip()))
                newfile.write(newline)
            else :
                newfile.write(s)
        print "Completed!"

    # Change the ID of residue, or index of residue
    def resSeqChanger(self, inline, resseq):
        resseqstring = " "*(4- len(str(resseq)))+str(resseq)
        newline = inline[:22] + resseqstring + inline[26:]
        return newline

    def atomSeqChanger(self, inline, atomseq) :
        atomseqstring = " " * (5 - len(str(atomseq))) + str(atomseq)
        newline = inline[:6] + atomseqstring + inline[11:]
        return newline

    def resNameChanger(self, inline, resname) :
        resnamestr = " " * ( 4 - len(str(resname) ) ) + str(resname)
        newline = inline[:17] + resnamestr + inline[20:]
        return newline

    def chainIDChanger(self, inline, chainid) :
        newline = inline[:21] + str(chainid) + inline[22:]
        return newline

class GenerateTop :
    def __init__(self):
        pass

    def gmxTopBuilder(self, frcmodFile, prepfile, PDBFile,outputName, ionName, ionNum, solveBox=None, boxEdge=12,
                      FField=["AMBER99SB", ], verbose=True):
        '''
        provide a coordination file, output a amber and gromacs topology file.

        required parameters

        :param frcmodFile: not required if the molecule is a protein, or DNA, RNA
        :param prepfile: not required if the molecule is a protein, or DNA, RNA
        :param PDBFile: input, coordination file
        :param outputName: output, output file name
        :param ionName: optional
        :param ionNum: optional
        :param solveBox: optional, type of solvation box, default is TIP3PBOX
        :param boxEdge: optional
        :param FField: a set of amber force fields, generally, we could choose
            99sb, 99sbildn, gaff, 12, 14sb
        :param verbose:
        :return: None
        '''

        # check AMBERHOME PATH
        AMBERHOME = sp.check_output("echo $AMBERHOME", shell=True)
        AMBERHOME = AMBERHOME[:-1] + "/"
        if verbose :
            print("Here is amber home path: " + AMBERHOME + " \n\n")

        # multiple FF supoorted here
        leapin = open("leapIn.in", 'w')
        for ff in FField:
            # create tleap input file
            if "gaff" in ff:
                leapin.write("source leaprc.gaff \n")
            elif "ildn" in ff or "ILDN" in ff:
                leapin.write("source oldff/leaprc.ff99SBildn  \n")
            elif ff == "AMBER99SB" or "amber99sb" == ff:
                leapin.write("source oldff/leaprc.ff99SB  \n")
            elif ff == "AMBER14SB" or "14SB" in ff:
                leapin.write("source ff14SB  \n")
            elif "leaprc." in ff:
                leapin.write("source %s  \n" % ff)
            else:
                print "Load Force Field File Error! \nExit Now!"
                sys.exit(1)

        # load amber frcmod and prep files
        if frcmodFile :
            leapin.write("loadamberparams  " + frcmodFile + "  \n")
        if prepfile :
            leapin.write("loadamberprep " + prepfile + "  \n")

        # prepare PDB file and load it
        if ".pdb" in PDBFile:
            leapin.write("pdb = loadPDB  " + PDBFile + "  \n")
        elif ".mol2" in PDBFile :
            # convert the mol2 file to pdb file using obabel
            job = sp.Popen("obabel %s -O %s"%(PDBFile, PDBFile[:-4]+"pdb"), shell=True)
            job.communicate()
            leapin.write("pdb = loadpdb " + PDBFile[:-4]+"pdb" + "  \n")
        elif len(PDBFile) >= 2:
            leapin.write("pdb = sequence{ ")
            for item in PDBFile:
                leapin.write(item + " ")
            leapin.write(" } \n")
        else:
            print "Loading PDB file or Sequence file error!"
            sys.exit(1)

        # save a amber off file of the molecule
        leapin.write("saveoff pdb %s.lib \n" % outputName)

        # add counter ions and solvate solute into water box
        if ionNum[0] > 0 and ionName[0] != "X+":
            if len(ionNum) == len(ionName):
                for i in range(len(ionNum)):
                    leapin.write("addions2 pdb %s %d \n" % (ionName[i], ionNum[i]))
            else:
                print "\nAdd ions not successful!\n"
        else:
            "\nNot adding any ions!\n"
        if solveBox :
            if boxEdge :
                leapin.write("solvatebox pdb %s %f \n" % (solveBox, boxEdge))
            else:
                print "\nBOX size not correctly set.\nExit Now!\n"
                sys.exit(1)
        else:
            print "\nNot setting simulation box!\n"

        # check object
        leapin.write("check pdb \n")
        leapin.write("savepdb pdb %s  \n" % (outputName + ".pdb"))
        leapin.write("saveoff pdb %s  \n"%(outputName+".lib"))
        leapin.write("saveamberparm pdb %s.prmtop %s.prmcrd \n" % (outputName, outputName))
        leapin.write("quit \n")
        leapin.close()

        if verbose :
            print("Generating a leap input file for tleap topology processing.")

        # run tleap
        out = sp.check_output("%s/bin/tleap -f leapIn.in  \n" % AMBERHOME, shell=True)
        if verbose :
            print(out)

        # convert AMBER format to GMX format
        time.sleep(1)
        out = sp.check_output("Amb2gmx.pl --prmtop %s.prmtop --crd %s.prmcrd --outname gmx_%s " \
                              % (outputName, outputName, outputName), shell=True)
        if verbose :
            print(out)
            print("GMX topology created")

    def runObabel(self, obabelexe, input, output, verbose=True):
        status = sp.check_output("%s %s -O %s "%(obabelexe, input, output))
        if verbose :
            print(status)
        return(1)

    def removeMolInfor(self, outputName, topFileName=None, verbose=True):
        # generate itp file from gmx.top file
        itpfile = open(outputName + ".itp", 'w')
        if not topFileName :
            topFileName = "gmx_" + outputName + ".top"

        with open(topFileName) as lines:
            start = 0
            stop = 0
            for s in lines:
                if "[ moleculetype" in s:
                    start = 1
                if "[ system" in s:
                    stop = 1

                if start and not stop:
                    if "solute" in s:
                        itpfile.write("%-6s            3 \n" % outputName)
                    else:
                        itpfile.write(s)
        itpfile.close()
        if verbose :
            print("ITP file for %s generated! " % outputName)
        return(1)

    def runAntechamber(self, infile, netCharge, antechamber="sample_antechamber.sh"):
        '''
        run antechamber to generate RESP am1/bcc atomic charges
        :param netcharge: input, number of netcharges of the molecule
        :param antechamber: input, optional, a sample antechamber shell script
        :return: the file name of the antechamber shell script
        '''
        if not os.path.exists(antechamber) :
            with open(antechamber, 'w') as tofile :
                content = '''
if [ $# -ne 2 ]; then
echo "Usage: input RED III.1  output, charge & spin & residuNAME read from Mol2"
echo "Format: file_in(Mol2); atom_type(gaff or amber)"
exit
fi

antechamber -i $1 -fi mol2 -o prep.$2 -fo prepi -at $2 -pf y -s 2
#antechamber -i $1 -fi mol2 -o prep.$2 -fo prepi -at $2 -pf y -s 2 -c resp
parmchk -i prep.$2 -f prepi -o frcmod.$2 -a Y
grep need frcmod.$2

if test -f leap.in; then
   rm leap.in
fi

echo "============================================="
echo "dmfff = loadamberparams frcmod.$2" >> leap.in
echo "loadamberprep prep.$2" >> leap.in

echo "prepared run_parm.sh for tleap"
echo "tleap -s -f ${AMBERHOME}dat/leap/cmd/leaprc.ff99 -f leap.in"
echo "tleap -s -f ${AMBERHOME}dat/leap/cmd/leaprc.gaff -f leap.in"
echo "general AMBER : leaprc.gaff"
echo "AMBER : leaprc.ff90"
                '''
                tofile.write(content)
        # run antechamber here, may first generate a sh script and then wrap it up
        if not netCharge :
            try:
                nc = SummaryPDB(pdbfile=infile)
                netcharge = nc.netCharges(inputMol=infile)
            except :
                print("Getting netcharge error!")
                netcharge = 0
        else :
            netcharge = 0

        tofile = open("antechamber.sh", 'w')
        with open(antechamber) as lines:
            for s in lines:
                if len(s.split()) > 0 and s.split()[0] == "antechamber":
                    tofile.write(
                        "antechamber -i $1 -fi mol2 -o prep.$2 -fo prepi -at $2 -pf y -s 2 -c bcc -nc %d \n" % netcharge)
                else:
                    tofile.write(s)
        return("antechamber.sh")

    def arguments(self):
        '''
        Parse arguments in the command line mode.

        :return: args, a parser.argument object
        '''

        # go current directory, pwd
        os.chdir(os.getcwd())

        if 0 :
            frcmodFile = "fit_R2.frcmod"
            PDBFile = "trpcage.pdb"
            outputName = "TRP"
            FField = "AMBER99SB"
            self.gmxTopBuilder(frcmodFile, PDBFile, outputName, FField)

        else:
            d = '''
            GMX-AMBER, topology builder.\n
            Provide a PDB file, or sometimes prep and frcmod files, itp and prmtop file will be given.
            Copyright @ Liangzhen Zheng, contact astrozheng@gmail.com for any technique support. \n
            Examples:
              python autoMD.py gentop -h
              python autoMD.py gentop -inp your.pdb -out yourpdb -ff AMBER14
              python autoMD.py gentop -aaseq ACE ALA CYS ALA HIS ASN NME -ff AMBER99SBildn -out peptide
              python autoMD.py gentop -inp lig.pdb -out lig -ff gaff -prep amber.prep.lig -frcmod frcmod.log
              python autoMD.py gentop -inp cplx.pdb -out cplx_box -ff gaff ildn
                                      -bt TIP3PBOX -d 1.2 -prep amber.prep.lig -frcmod frcmod.log
              python autoMD.py gentop -inp ligand.mol2 -resp True -nt 0 -out LIG
            '''
            parser = argparse.ArgumentParser(description=d, formatter_class=RawTextHelpFormatter)
            parser.add_argument("-inp", type=str, default=None,
                                help="The PDB file or a mol2 file for topology generating. Default is None\n")
            parser.add_argument("-prep", type=str, default="amberff.prep",
                                help="Prep file (amber format) stored coordinates and charges.\n"
                                     "Default is None.")
            parser.add_argument("-frcmod", type=str, default="frcmod",
                                help="The additional parameters, stored in frcmod file (Amber format).\n")
            parser.add_argument("-out", type=str, default="OUT", help="The output name. Default is OUT.\n")
            parser.add_argument("-ion", type=str, nargs="+", default=["X+", ],
                                help="The ions to be added to the system. Mutiple Ions supported. \n"
                                     "Options are Na+, Cl-\n"
                                     "Default is X+. \n")
            parser.add_argument("-nion", type=int, nargs="+", default=[-1, ],
                                help="Number of ions to be added. Default is -1. \n"
                                     "Number of args must be the same as -ion. \n")
            parser.add_argument("-bt", type=str, default=None,
                                help="Box type, ff you wish to solvate the mol in water box, you should\n"
                                     "provide an option: NA, TIP3PBOX, TIP4PEW, or TIP5PBOX. \n"
                                     "Default choice is TIP3PBOX. \n")
            parser.add_argument("-size", type=float, default=0,
                                help="The size of your solvation box. Default is 0 angstrom. \n")
            parser.add_argument("-ff", type=str, nargs='+', default=["AMBER99SB", ],
                                help="The force field for simulation. Multiple force filed files were supported.\n"
                                     "Options including: 99SB, 99SBildn, AMBER14, gaff\n"
                                     "or use leaprc.xxx files here. Default is AMBER99SB. \n")
            parser.add_argument("-aaseq", type=str, nargs='+', default=None,
                                help="Amino acid sequences in case no PDB file provide.\n"
                                     "Default is None. It is an optional argument.\n")
            parser.add_argument("-resp", default=False, type=bool,
                                help="If it is a small molecule, whose atomic charges \n"
                                     "not defined, neither the bonding angles, we could \n"
                                     "use RESP with bcc AM1 to determine charges and \n"
                                     "prepare parameters. Options: True, False. \n"
                                     "Defualt is False. \n")
            parser.add_argument("-nc", default=None, type=int,
                                help="Net charge. This argument only works with \n"
                                     "the -resp argument. Default is 0. \n")

            # args include the options in the argparser
            args, unknown = parser.parse_known_args()
            if len(sys.argv) < 3 :
                parser.print_help()
                print("Number of arguments are not correct! Exit Now!")
                sys.exit(1)

            return(args)

class ExtractPDB :
    def __init__(self):
        pass

    def extract_pdb(self,filename, structname, first_n_frame):
        '''
        from a big pdb file to extract single PDB structure file

        :param filename:
        :param structname:
        :param first_n_frame:
        :return:
        '''

        lines = open(filename)
        file_no = 1
        pdbline = open(structname+'_'+str(file_no)+'.pdb','w')

        for s in lines :

            if  "MODEL" in s :
                if file_no != 1 :
                    pdbline = open(structname+'_'+str(file_no)+'.pdb','w')
                pdbline.write(s)
                print "Start Model " + str(file_no)
                file_no += 1

            elif "ATOM" == s.split()[0] or "TER" == s.split()[0] :
                pdbline.write(s)

            elif "ENDMDL" in s :
                pdbline.write(s)
                pdbline.close()

            else :
                pass

            if first_n_frame != 0 and file_no > first_n_frame + 1 :
                print "Extract only the first "+str(first_n_frame)
                break

        print "Finished!\n\n"

    def extract_frame(self, filename, structname, no_frames=[]):
        lines = open(filename)
        print "The models of the pdb file is : "
        for s in lines :
            if "MODEL" in s :
                print "    "+s[:-1]
        lines.close()

        if len(no_frames) :
            try :
                print "Which frames would you want to extract ? "
                frames = raw_input("Input the frame number(s) here (multi numbers are accepted):  ")
                frame_list = frames.split()
            except :
                print("You haven't select correct frames.")
        else :
            frame_list = no_frames

        for frame_no in frame_list :

            lines  = open(filename)
            condition = False
            for s in lines :
                if "MODEL" in s and int(frame_no) == int(s.split()[1])  :
                    newline = open(structname+"_"+str(frame_no)+".pdb","w")
                    newline.write(s)
                    condition = True
                elif "ATOM" in s and condition :
                    newline.write(s)
                elif condition and "ENDMDL" in s :
                    condition = False
                elif "MODEL" in s and int(frame_no)+1 == int(s.split()[1]) :
                    condition = False
                    break
                else :
                    condition = False
            newline.close()
            lines.close()

        print "Finished writing frames to separated files!\n\n"
        return(1)

    def printinfor(self):
        print "What would you like to do now ?"
        print "  1. Extract all the frames from the input pdb file;"
        print "  2. Extract selected frames from the input file;"
        print "  3. Extract the first N frames from the input file;"
        print "  4. Do nothing and exit now. "

    def indexCoord(self, filename):
        '''
        provide a large coordination file, eg, a pdb file or a mol2 file,
        return the indexing of the frames

        :param filename: the file name of a large multiple-frames coordination file
            either a pdb, a pdbqt or a mol2 file
        :return: the indexing of the first line of a multiple-frame file
        '''
        indexing = []
        lineNumber=-1

        extention = filename.split(".")[-1]

        with open(filename) as lines :
            if extention in ['pdb', 'pdbqt'] :
                for s in lines :
                    lineNumber += 1
                    if len(s.split()) > 1 and "MODEL" == s.split()[0] :
                        indexing.append(lineNumber)
            elif extention in ['mol2'] :
                for s in lines :
                    lineNumber += 1
                    if "@<TRIPOS>MOLECULE" in s :
                        indexing.append(lineNumber)
            else:
                print("Only a pdb file, a pdbqt file or a mol2 file supported.")
                sys.exit(0)

        return(indexing)

    def extract_all(self, filename, structname):
        '''
        extract all the frames into seperated mol2 (or, pdb and pdbqt) files

        :param filename: the multiple-frames mol2 file
        :param structname: the prefix of the extracted separated frames
        :return: None
        '''
        extension = filename.split('.')[-1]
        if extension in ['pdb', 'pdbqt', 'mol2'] :
            try :
                # try to loop over the file to count number of lines in the file
                totalLineNum = sum(1 for line in open(filename))
            except IOError :
                totalLineNum = 0

            # if in the file filename, the file is not empty,
            # start to extract frames
            if totalLineNum :
                structFirstLineIndex = self.indexCoord(filename)
                # at the end of file, provide a sudo-next frame start line index
                structFirstLineIndex.append(totalLineNum)

                for i in range(len(structFirstLineIndex))[1:] :
                    start_end = [structFirstLineIndex[i-1], structFirstLineIndex[i]]

                    with open(structname+"_"+str(i)+"."+extension, 'wb') as tofile :
                        # extract the specific lines from the large multiple-frame
                        # file to write to a new file.
                        for lndx in range(start_end[0]+1, start_end[1]+1) :
                            tofile.write(linecache.getline(filename, lndx))
            else :
                print("File %s is empty. Couldnot extract frames. " % filename)
        else :
            print("PDB, PDBQT, or MOL2 file is required. ")
            sys.exit(0)

        print("Extracting all frames in mol2 file completed. ")

    def arguments(self):
        # PWD, change directory to PWD
        os.chdir(os.getcwd())

        d = '''
        This script try to extract PDB frames from a long trajectory file from GMX trjconv or g_cluster.
        Any questions, contact Liangzhen Zheng, astrozheng@gmail.com

        Examples :
        Get help information
        python autoMD.py extract -h
        Extract frames in a multiple-frames PDB file, which contains "MODEL 1", "MODEL 2"

        PDB_Extractor.py extract -i md.pdb -o splitted

        after which interactive mode is entered and you are promoted to choose one of the 4 options:
            1. Extract all the frames from the input pdb file;
            2. Extract selected frames from the input file;
            3. Extract the first N frames from the input file;
            4. Do nothing and exit now.
        '''
        parser = argparse.ArgumentParser(description=d, formatter_class=RawTextHelpFormatter)
        parser.add_argument('-i', '--input', type=str, default="yourpdb.pdb",
                            help="Input PDB file name. Default is yourpdb.pdb.")
        parser.add_argument('-o', '--output', type=str, default="out",
                            help="Output file format. Default is out_*.pdb.")
        parser.add_argument('-m', '--allMol2', type=bool, default=False,
                            help='Extract all the frames in a mol2 file. \n'
                                 'Options: True, False \n'
                                 'Default is False.')
        options, args = parser.parse_known_args()

        if len(sys.argv) < 2:
            parser.print_help()
            sys.exit(1)
        else:
            parser.print_help()

        return(options)

class SummaryPDB :
    def __init__(self, pdbfile, aminoLibFile=''):
        self.pdbfile = pdbfile

        resShortName = {}
        try :
            with open(aminoLibFile) as lines:
                for s in lines:
                    if "#" not in s:
                        resShortName[s.split()[2]] = s.split()[3]
            self.resShortName = resShortName
        except IOError :
            self.resShortName = {}

    def netCharges(self, inputMol, ligName=None):
        '''
        netCharges determine the net charge of a molecule, given a pdb file.
        if the input molecule is a pdbqt (Vina input), or pqr (APBS input type),
           a molecule name (residue code) is need.
        if a mol2 file provide, only the @<TRIPOS>ATOM field will be used, no ligand name required.

        other formats are not supported.

        last column of a pdbqt, pqr and mol2 file generally will be the atomic charge field,
          otherwise, a ValueError exception will be rasied.

        :param inputMol: input file with atomic charges in the last column
        :return int netCharge
        '''
        extension = inputMol.split(".")[-1]
        netCharge = 0.0

        with open(inputMol) as lines :
            if extension in ['pdbqt', 'pqr']:
                for s in lines :

                    if s.split()[0] in ["ATOM", "HETATM"] and len(s.split()) > 5 :
                        try :
                            if ligName and ligName in s :
                                netCharge += float(s.split()[-1])
                            elif not ligName :
                                netCharge += float(s.split()[-1])
                        except ValueError :
                            netCharge += 0.0
                            print("Last column in %s is not a float point charge value."%inputMol)
            elif extension in ['mol2'] :
                condition = False
                for s in lines :
                    if len(s.split()) and "@" in s :
                        if "<TRIPOS>ATOM" in s :
                            condition = True
                        else :
                            condition = False
                    elif condition and len(s.split() ):
                        try:
                            netCharge += float(s.split()[-1])
                        except ValueError :
                            netCharge += 0.0
                            print("Last column in %s is not a float point charge value." % inputMol)
            else :
                print("The file extention is not recognized. "
                      "\nPlease provide a pdbqt, pqr, or a mol2 file.")

        return(int(netCharge))

    def centerOfMass(self, inputMol, atomNdx, obabelexe='obabel', molBox=False):
        '''
        Given a file (preferable PDB file format),
            if not, the file will be convert into a pdb file,
        and selected atom sequence number,
        determine the COM of the coordinates

        :param inputMol: a file, the input coordinates
        :param atomNdx: a list, atom sequence number
        :return: two lists, center of mass, format [0.0, 0.0, 0.0]
                            boxsize, format [10.0, 10.0, 10.0]
        '''
        pdb = inputMol
        if inputMol.split(".")[-1] not in ['pdb', 'pdbqt'] :
            spdb = GenerateTop()
            spdb.runObabel(obabelexe, inputMol, inputMol+".pdb")
            pdb = inputMol + ".pdb"
        coordinates = []

        with open(pdb) as lines :
            for s in lines :
                if "ATOM" == s.split()[0] or "HETATM" == s.split()[0] and s.split()[1] in atomNdx :
                    crd_list = []
                    crd_list.append(float(s[30:38].strip()))
                    crd_list.append(float(s[38:46].strip()))
                    crd_list.append(float(s[46:54].strip()))
                    coordinates.append(crd_list)

        coordinates = np.asarray(coordinates)

        xcenter = np.mean(coordinates[:, 0])
        ycenter = np.mean(coordinates[:, 1])
        zcenter = np.mean(coordinates[:, 2])

        if molBox :
            xsize = 2 * max(np.max(coordinates[:, 0]) - xcenter,
                        np.abs(np.min(coordinates[:,0])-xcenter))
            ysize = 2 * max(np.max(coordinates[:, 1]) - ycenter,
                        np.abs(np.min(coordinates[:, 1]) - ycenter))
            zsize = 2 * max(np.max(coordinates[:, 2]) - zcenter,
                        np.abs(np.min(coordinates[:, 2]) - zcenter))
        else :
            xsize, ysize, zsize = 100, 100, 100
        return([xcenter, ycenter, zcenter], [xsize, ysize, zsize])

    def details(self, verbose=False):
        chains = []
        resNdx = defaultdict(list)
        resName= defaultdict(list)
        resAtom= defaultdict(list)
        resNameNdx = {}
        with open(self.pdbfile) as lines :
            for s in lines :
                if len(s.split()) > 0 and s.split()[0] in ["ATOM","HETATM"] and "TER" not in s :
                    if s[21] not in chains:
                        chains.append(s[21])

                    ## residue index information
                    if s[21] not in resNdx.keys() :
                        resNdx[s[21]] = []
                    if s[22:26].strip() not in resNdx[s[21]] :
                        resNdx[s[21]].append((s[22:26].strip()))

                    ## residue name information
                    if s[21] not in resName.keys():
                       resName[s[21]] = []
                    if s[17:20].strip() + "_" + s[22:26].strip() not in resName[s[21]] :
                       resName[s[21]].append(s[17:20].strip() + "_" + s[22:26].strip())

                    # pdb file atoms name in each residue
                    resId = (s[22:26].strip() + '_' + s[21]).strip()
                    if resId not in resAtom.keys() :
                        resAtom[resId] = []
                    resAtom[resId].append(s[12:16].strip())

                    # residue index and name hash map
                    resNameNdx[s[22:26].strip()+'_' + s[21].strip()] = s[17:20].strip()

        ## print some information
        if verbose :
            print "\nNumber of chains in this pdb file : ", len(chains), ". They are ", chains

            print "\nNumber of residues in each chain: "
            for chain in chains :
                print "Chain ", chain, " ",len(resNdx[chain])

            print "\nResidue Names for each chain are : "
            for chain in chains :
                print "For chain ", chain
                for i in range(10) :
                    print resNdx[chain][i], "  ", resName[chain][i]
                print "......"
                for j in range(10) :
                    print resNdx[chain][-10+j], "  ", resName[chain][-10+j]
        return chains, resNdx, resName, resAtom, resNameNdx

    def summary(self, chain, verbose=False):
        chains, resNdx, resName, resAtom, resNameNdx = self.details()

        proteinSequence = {}
        noProteinResNdx = defaultdict(list)
        noProteinResName= defaultdict(list)
        missingResNdx = defaultdict(list)
        fullResNdx = defaultdict(list)
        #for chain in chains :
        #resNdxNameDict = {}
        if len(resNdx[chain]) != len(resName[chain]) :
            print "Error in function SummaryPDB.details. " \
                  "\nNumber of index is different with number of residues."
            print "Exit Now!"
            sys.exit(1)

        proResNdx = []
        for i in range(len(resNdx[chain])) :
            ## get index for all the protein residues in a specific chain
            if resName[chain][i].split("_")[0] in self.resShortName.keys() :
                proResNdx.append(int(resNdx[chain][i]))

                #resNdxNameDict[int(resNdx[chain][i])] = resName[chain][i].split("_")[0]
            else :
                ## non-protein residues information
                if chain not in noProteinResName.keys() :
                    noProteinResName[chain] = []
                if chain not in noProteinResNdx.keys() :
                    noProteinResNdx[chain]  = []

                noProteinResName[chain].append(resName[chain][i])
                noProteinResNdx[chain].append(resNdx[chain][i])

        ## get protein chain sequence
        startNdx = proResNdx[0]
        finalNdx = proResNdx[-1]
        #print startNdx,finalNdx, chain, " chain starting and end index"

        fullNdx = range(startNdx, finalNdx+1)
        fullResNdx[chain] = fullNdx
        #print fullNdx
        for i in range(len(fullNdx) ):
            residueid = str(fullNdx[i]) + "_" + chain
            #print residueid
            #print i, fullNdx[i]
            if chain not in proteinSequence.keys() :
                proteinSequence[chain] = ''

            if residueid in resNameNdx.keys() :
                if resNameNdx[residueid] not in self.resShortName.keys() :
                    proteinSequence[chain] += "_" + resNameNdx[residueid]
                else :
                    proteinSequence[chain] += self.resShortName[resNameNdx[residueid]]
            else :
                proteinSequence[chain] += "-"

                ## find missing residues' index
                if chain not in missingResNdx.keys() :
                    missingResNdx[chain] = []
                missingResNdx[chain].append(str(fullNdx[i]))
        #print chain, proteinSequence[chain]

        ## print some information
        if verbose :
            print "\nFull sequence of protein chains are: "
            #for chain in chains :
            print chain, "  ", proteinSequence[chain]

            print "\nSome missing protein residues in each chain "
            #for chain in chains :
            print chain, "  ", missingResNdx[chain]

            print "\nThe information of the non-protein residues here "
            #for chain in chains :
            print chain, "  ", noProteinResName[chain]

            print "\nThe sequence of the full protein chain is: "
            print "Chain  ", chain
            sections = math.ceil(float(len(proteinSequence[chain])) / 20.0)

            for i in range(int(sections - 1)):
                print fullResNdx[chain][i * 20], " " * (-len(str(fullResNdx[chain][i * 20])) + 11), fullResNdx[chain][
                    i * 20 + 10]
                print proteinSequence[chain][i * 20: i * 20 + 10], " ", proteinSequence[chain][i * 20 + 10: i * 20 + 20]
            print fullResNdx[chain][i * 20 + 20]
            print proteinSequence[chain][i * 20 + 20:]

        return proteinSequence, missingResNdx, noProteinResNdx, noProteinResName, fullResNdx

    def missingRes(self, chain, fastaSeq, verbose=False):
        ## find the missing residues sequences in the pdb file
        #chains, resNdx, resName, resAtom, resNameNdx = self.details()
        proteinSequence, missingResNdx, noProteinResNdx, noProteinResName, fullResNdx = self.summary(chain)
        #print fullResNdx[chain]
        trueResName = defaultdict(list)
        ## full fasta sequences here :
        fastaseq = ''
        if os.path.isfile(fastaSeq):
            with open(fastaSeq) as lines:
                for s in lines:
                    if '>' not in s :
                        #print s #strip()
                        fastaseq += s.strip()
        #print fastaseq
        else :
            fastaseq = fastaSeq

        # print protein missing residue information
        startndx = fullResNdx[chain][0]
        for index in missingResNdx[chain] :
            if chain not in trueResName.keys() :
                trueResName[chain] = []
            trueResName[chain].append(fastaseq[int(index)-startndx]+'_'+index)

        ranges = []
        siteRange = []
        missedSeqsList = defaultdict(list)
        missedRSeq = ''
        brokenSeq = proteinSequence[chain]
        for i in range(len(brokenSeq)) :
            currentResNdx = i + int(startndx)
            if brokenSeq[i] != '-' and len(siteRange) == 0 :
                pass
            elif brokenSeq[i] == '-' :
                siteRange.append(currentResNdx)
                missedRSeq += fastaseq[i]
                #print missedRSeq
            elif brokenSeq[i] != '-' and len(siteRange) > 0 :
                ranges.append([siteRange[0],siteRange[-1]])
                #missedRSeq.append(seq[i])

                missedSeqsList[str(siteRange[0])+"_"+str(siteRange[-1])] = missedRSeq
                missedRSeq = ''
                siteRange = []
            else :
                pass
        #print missedSeqsList

        if verbose :
            ## print full sequence information
            print "The fasta sequence of the full protein chain is: "
            print "Chain  ", chain
            sections = math.ceil(float(len(fastaseq)) / 20.0)

            for i in range(int(sections-1)) :
                print fullResNdx[chain][i*20], " "*(-len(str(fullResNdx[chain][i*20]))+11), fullResNdx[chain][i*20+10]
                print fastaseq[i*20 : i*20+10], " ",fastaseq[i*20+10 : i*20+20]
            print fullResNdx[chain][i*20+20]
            print fastaseq[i*20+20:]

            ## print the true residue name
            print "\nThe missing residues are :"
            print chain,"  ", trueResName[chain]

            ## information about the missing sequences
            for ndxrange in missedSeqsList.keys() :
                print (("%12s   %-s ")%(ndxrange, missedSeqsList[ndxrange]))


        return fastaseq, missedSeqsList, fullResNdx

    def extendMissRange(self, chain, fastaSeq, extend=10, verbose=False):
        extMissSeq = {}
        fastaseq, missedSeqs, fullResNdx = self.missingRes(chain,fastaSeq)
        startndx = int(fullResNdx[chain][0])
        for ndxrange in missedSeqs.keys() :
            fixSeq = fastaseq[(int(ndxrange.split("_")[0])-extend-startndx):(int(ndxrange.split("_")[1]) + extend-startndx+1)]
            newrange = str((int(ndxrange.split("_")[0])-extend)) + "_" +\
                str(int(ndxrange.split("_")[1]) + extend)
            extMissSeq[newrange] = fixSeq

        if verbose :
            print "\nExtend %d residues at broken part for chain %s : "%(extend,chain)
            for newrange in extMissSeq.keys() :
                print(("%12s   %-s")%(newrange, extMissSeq[newrange]))

        return extMissSeq

class FixPBD :
    def __init__(self):
        pass

    def pdbDownloader(self, pdbCode, pdbout):
        if not os.path.exists(pdbCode + '.pdb'):
            try :
                source = urllib2.urlopen("http://www.rcsb.org/pdb/files/" + pdbCode + ".pdb")
                with open(pdbCode + '.pdb', 'wb') as target:
                    target.write(source.read())
            except urllib2.URLError, e :
                print e.args

        return(1)

    def removeRegions(self, filename, residuesNdx, chain, pdbout="temp_1.pdb", verbose=False):
        '''
        input a pdbfile, remove the selected residues

        :param filename: input pdbfile
        :param residuesNdx: the range of residue index for deleting
        :param chain: which chain to perform delete
        :param pdbout: the output pdb file name
        :return: 1
        '''
        tofile = open(pdbout,'w')
        with open(filename) as lines :
            for s in lines :
                if s.split()[0] in ["ATOM", "HETATM" ] and s[21] == chain and int(s[22:26].strip()) in residuesNdx :
                    pass
                else :
                    tofile.write(s)
        tofile.close()

        return(1)

    def addModeledRegions(self, basepbd, modeledpdb,
                          modelNdx, removeNdx, chain,
                          pdbout="temp.pdb",
                          verbose=False
                          ):
        tofile = open(pdbout, 'w')
        with open(basepbd) as lines:
            for s in lines :
                if s.split()[0] in ["ATOM", "HETATM" ] and s[21] == chain and int(s[22:26].strip()) < removeNdx[0] :
                    tofile.write(s)
                else :
                    pass

        ## addd modeled pdb here
        with open(modeledpdb) as lines :
            for s in lines :
                if s.split()[0] in ["ATOM", "HETATM"] and int(s[22:26].strip()) in modelNdx :
                    tofile.write(s)
                else :
                    pass

        ## add the following part of the original PDB
        with open(basepbd) as lines:
            for s in lines :
                if s.split()[0] in ["ATOM", "HETATM"] and s[21] == chain and int(s[22:26].strip()) > removeNdx[-1]:
                    tofile.write(s)
                else :
                    pass

        tofile.close()
        return(1)

    def addMissingAtoms(self, pdbin, pdb2pqrPy,
                        ff='amber', pqrout='output.pqr',
                        addhydrogen=True,
                        ph=7.0, verbose=False ):
        '''
        add missing heavy atoms, and assign hydrogens according to the ph

        :param pdbin:
        :param pdb2pqrPy:
        :param ff:
        :param pqrout:
        :param addhydrogen:
        :param ph:
        :return:
        '''
        if addhydrogen :
            cmd = "python %s -ff %s --with-ph %f --ffout %s --chain %s %s " %\
                  (pdb2pqrPy, ff, ph, ff+"_"+pqrout, pdbin, pqrout)
        else :
            cmd = "python %s -ff %s --ffout %s --chain %s %s " %\
                  (pdb2pqrPy, ff, ff+"_"+pqrout, pdbin, pqrout)
        job = sp.Popen(cmd, shell=True)
        job.communicate()
        job.kill()

        return(ff+"_"+pqrout)

    def addLoopsAutoModeller(self, pdbCode, chain1,
                         alignCode, chain2,
                         loopRange,
                         verbose=False):
        '''
        model missing loop region
        :param pdbCode:
        :param chain1:
        :param alignCode:
        :param chain2:
        :param loopRange:
        :param verbose: print more information
        :return: the final model pdb file
        '''
        if MODELLER_EXIST :
            if verbose :
                print("MODELLER exist, perform alignment and loop refinement.")
                log.verbose()

            for pdb in [pdbCode]+[alignCode] :
                if not os.path.exists(pdb+'.pdb') :
                    self.pdbDownloader(pdb, pdb+'.pdb')

            env = environ()
            # directories for input atom files
            env.io.atom_files_directory = ['.', '../atom_files']

            aln = alignment(env)
            mdl_1 = model(env, file=pdbCode, model_segment=('FIRST:%s'%chain1, 'LAST:%s'%chain1))
            aln.append_model(mdl_1, align_codes=pdbCode, atom_files=pdbCode+'.pdb')
            mdl_2 = model(env, file=alignCode, model_segment=('FIRST:%s'%chain2, 'LAST:%s'%chain2))
            aln.append_model(mdl_2, align_codes=alignCode, atom_files=alignCode+'.pdb')
            aln.align2d()
            aln.write(file='alignment.ali', alignment_format='PIR')

            """ select missing loop residues for modeling only"""
            class MyModel(automodel):
                # overide select_atom function, select only some residues
                def select_atoms(self):
                    # Select residues loopRange[0] to loopRange[1] (PDB numbering)
                    return selection(self.residue_range(str(loopRange[0])+":"+str(chain1),
                                                        str(loopRange[1])+":"+str(chain1)
                                                        )
                                     )

            a = MyModel(env,
                        alnfile='alignment.ali',
                        knowns=alignCode,
                        sequence=pdbCode,
                        )
            # a.auto_align()
            # get an automatic loop building
            a.make()

            # obtain successfully modeled pdb file names
            ok_models = [x for x in a.outputs if x['failure'] is None]

            # Rank the models by DOPE score
            key = 'DOPE score'
            if sys.version_info[:2] == (2, 3):
                # Python 2.3's sort doesn't have a 'key' argument
                ok_models.sort(lambda a, b: cmp(a[key], b[key]))
            else:
                ok_models.sort(key=lambda a: a[key])

            # Get top model
            finalPDB = ok_models[0]
            if verbose :
                print("Top model: %s (DOPE score %.3f)" % (finalPDB['name'], finalPDB[key]))

            return(finalPDB)
        else :
            if verbose :
                print("No model built!")
            return(pdbCode+'.pdb')

class CleanPDB :
    '''
    Prepare the pdb file, add missing atoms, add hydrogens, add missing loops
    '''
    def __init__(self, filename, obabel='obabel'):
        self.pdb = filename
        # the input file is not a pdb, or pdbqt file, use obabel convert it
        if filename.split(".")[-1] not in ["pdb","pdbqt"] :
            gpdb = GenerateTop()
            gpdb.runObabel(obabelexe=obabel, input=filename, output=filename+".pdb")
            self.pdb = filename+".pdb"
    def extractFrame(self, frame=1, headersave=True) :
        '''
        extract a frame (default is first) of the pdb file

        :param frame:
        :param headersave:
        :return: the specific frame from a large pdb file
        '''
        extractPDB = ExtractPDB()
        extractPDB.extract_frame(self.pdb, self.pdb[:-5]+"_%d.pdb"%frame, no_frames=[frame])
        return(self.pdb[:-5]+"_%d.pdb"%frame)

    def removeHETATM(self, filename, hetatmsave=['WAT','HOH'], dropWater=True,
                     cleanedPDB='cleaned_', headersave=True,
                     selectedChains = [],
                     ):
        '''
        :param filename:
        :param hetatmsave:
        :param dropWater: remove water molecules
        :param cleanedPDB:
        :param headersave: save the pdb header information
        :param selectedChains:
        :return: the pdbfile after removing unnecessary information
        '''
        tofile = open(cleanedPDB+filename, 'wb')
        with open(filename) as lines :
            for s in lines :
                if len(s.split()) > 0 and \
                                s.split()[0] in ['ATOM', 'HETATM', 'TER', 'END', 'ENDMDL'] :
                    if dropWater :
                        if s[21] in selectedChains and s[17:20].strip() not in ['HOH', 'WAT'] \
                                and s[17:20].strip() in hetatmsave :
                            tofile.write(s)
                    else :
                        if s[21] in selectedChains and s[17:20].strip() in hetatmsave :
                            tofile.write(s)
                else :
                    if headersave :
                        tofile.write(s)

        return(cleanedPDB+filename)

class MolDocking :
    def __init__(self, method):
        self.method = method

    def runVina(self, vinaexe, receptor, ligand,
                output='result.pdbqt', logfile='result.log',
                ncpu=1, exhaustiveness=32,
                center=[0,0,0], sizes=[40,40,40],
                no_modes = 20, en_range = 5, seed=-1,
                ):
        try :
            job = sp.Popen('%s --receptor %s --ligand %s '
                           '--center_x %f --center_y %f --center_z %f '
                           '--size_x %f --size_y %f --size_z %f '
                           '--log %s --out %s --cpu %d '
                           '--exhaustiveness %d --num_modes %d '
                           '--energy_range %d --seed %d ' %
                           (
                               vinaexe, receptor, ligand,
                               center[0], center[1], center[2],
                               sizes[0], sizes[1], sizes[2],
                               logfile, output, ncpu,
                               exhaustiveness, no_modes,
                               en_range, seed),
                           shell= True
                           )
            job.communicate()
            job.terminate()
        except :
            print("Docking molecule %s to %s using %s failed. \n"
                  "Check your input and logfile.")

        return(1)

class PrepScripts :
    def __init__(self, qsubScriptSample):
        self.qsubMD = qsubScriptSample

    def prepareQsub(self,  outSH, cmds):
        tofile = open(outSH, 'wb')

        with open(self.qsubMD) as lines:
            for s in lines:
                if "#define_commends" in s:
                    for cmd in cmds :
                        tofile.write(cmd+" \n")
                else:
                    tofile.write(s)
        tofile.close()

    def RewriteAntechamber(self, nc, antechamber="sample_antechamber.sh"):
        tofile = open("antechamber.sh", 'wb')
        with open(antechamber) as lines:
            for s in lines:
                if len(s.split()) > 0 and s.split()[0] == "antechamber":
                    tofile.write(
                        "antechamber -i $1 -fi mol2 -o amberff.prep.$2 -fo prepi -at $2 -pf y -s 2 -c bcc -nc %d \n" % nc)
                else:
                    tofile.write(s)
        return 1

class AutoRunMD :
    def __init__(self, topFile, taskName, grompp, mdrun, verbose=True, qsub=False):
        self.top = topFile
        self.task= taskName
        self.verbose = verbose
        self.grompp = grompp
        self.mdrun = mdrun

    def modifyMDP(self, inputMDPFile, outputMDPFile, parameters):
        tofile = open(outputMDPFile, 'wb')
        with open(inputMDPFile) as lines :
            for s in lines :
                if len(s.split()) > 0 and s[0] != ";" :
                    if s.split()[0] in parameters.keys() \
                            and len(parameters[s.split()[0]]) > 0 :
                        tofile.write("%s    = ")
                        for item in parameters[s.split()[0]] :
                            tofile.write("%s "%str(item))
                        tofile.write(" \n")
                    else:
                        tofile.write(s)
                else :
                    tofile.write(s)
        return(outputMDPFile)

    def energyMinimization(self, emMDPFile, groFile, outTprFile, np=4, qsub=False):
        if self.verbose :
            print("Start energy minimization for task %s."%self.task)

        # generate GMX tpr file
        job = sp.check_output("%s -f %s -c %s -o %s -p %s "%
                              (self.grompp, emMDPFile, groFile, outTprFile, self.top),
                              shell=True )
        if self.verbose :
            print(job)
            print("Generating TPR file for EM completed.")

        # run MD with qsub or on the terminal
        cmd = "mpirun -np %d %s -deffnm %s " % (np, self.mdrun, outTprFile[-4:])
        if qsub :
            pscript = PrepScripts('sample_qsub.sh')
            pscript.prepareQsub('em_qsub.sh', cmd)
            job = sp.Popen("qsub em_qsub.sh", shell=True)
            job.communicate()

        else :
            job = sp.Popen(cmd, shell=True)
            job.communicate()

        return(outTprFile[-4:]+".gro")

    def runMD(self, MDPFile, groFile, outTprFile, mdpOptions,
              sampleScript='sample_qsub.sh',
              qsub=False, np=4,
              index="index.ndx",
              ):

        # prepare mdp file for MD
        current_mdp = MDPFile
        if len(mdpOptions.keys()) :
            self.modifyMDP(inputMDPFile=MDPFile,
                           outputMDPFile="modified_"+MDPFile,
                           parameters=mdpOptions)
            current_mdp = "modified_"+MDPFile

        # prepare tpr file for md
        job = sp.check_output("%s -f %s -c %s -o %s -p %s -n %s"
                              %(self.grompp, current_mdp, groFile, outTprFile, self.top, index),
                              shell=True)
        if self.verbose :
            print(job)

        # run md
        cmd = "mpirun -np %d %s -deffnm %s " % (np, self.mdrun, outTprFile[-4:])
        if qsub :
            pscript = PrepScripts(sampleScript)
            oscript = str(len(glob("./*_"+sampleScript)))+"_"+sampleScript
            pscript.prepareQsub(oscript, cmd)
            job = sp.Popen("qsub %s"%oscript, shell=True)
            job.communicate()
        else :
            job = sp.Popen(cmd+" &", shell=True)
            job.communicate()

        return(outTprFile[-4:]+".gro")

def runExtract() :
    '''
    run ExtractPDB class
    extract frames from a multiple-frames pdb file or mol2 file
    :return: None
    '''
    pdb = ExtractPDB()
    args = pdb.arguments()

    pdbfile = args.input
    structname = args.output

    if args.allMol2 :
        # if we choose to extract all frames in a mol2 file
        pdb.extract_all(args.input, args.output)
    else :
        pdb.printinfor()
        command = raw_input("Your choice:  ")
        while command in ['0', '1', '2', '3', '4']:
            if command == "1":
                pdb.extract_pdb(pdbfile, structname, 0)
                command = '0'
            elif command == "2":
                pdb.extract_frame(pdbfile, structname)
                command = '0'
            elif command == "3":
                fnframe = raw_input("  The first N frames output : N = ")
                pdb.extract_pdb(pdbfile, structname, int(fnframe))
                command = '0'
            elif command == "4" or command not in "01234":
                print "\nExit now. \nHave a nice day!"
                command = '5'
                sys.exit(0)
                # return None
            elif command == "0":
                pdb.printinfor()
                command = raw_input("Your choice:  ")

def runGenTop() :
    '''
    instant a gmxTopBuilder class, and input some parameters,
        to generate amber and gromacs topology and coordinates

    if a small ligand provide, if not prep and frcmod exist,
        then am1/bcc charge based on amber antechamber sqm method
        will be calculated. in this process, tleap will be used to
        get bonded and non-bonded parameters
    :parameter: None
    :return: None
    '''
    top = GenerateTop()
    args = top.arguments()

    try :
        os.environ["AMBERHOME"] = "/home/liangzhen/software/amber14/"
    except :
        os.system("export AMBERHOME=/home/liangzhen/software/amber14/")

    # define the input coordinates information, a pdb, or mol2, or a sequence
    # of amino acids
    if args.aaseq:
        structure = args.aaseq
    elif args.inp:
        structure = args.inp
    else:
        structure = ''
        print('Error: No input structure defined.\n Exit Now!')
        sys.exit(0)

    if not args.resp:
        # not calculate atomic charges
        top.gmxTopBuilder(args.frcmod, args.prep,
                          structure, args.out,
                          args.ion, args.nion,
                          args.bt, args.size,
                          args.ff)
    else:
        # prepare resp charges for small ligand with antechamber
        # give a netcharge: args.nc
        sh = top.runAntechamber(args.nc)

        try:
            # set env for AMBERHOME
            os.environ["AMBERHOME"] = "/home/liangzhen/software/amber14/"
            job = sp.Popen("sh %s %s %s" % (sh, args.inp, args.out), shell=True)
            job.communicate()
        except :
            print("Antechamber generateing RESP and atomic parameters error! \n "
                  "Double check your input structure.")
            sys.exit(1)

        top.gmxTopBuilder("frcmod." + args.out,
                          "prep." + args.out,
                          structure, args.out,
                          args.ion, args.nion,
                          args.bt, args.size,
                          args.ff)

    if args.ion[0] == "X+" or args.bt == "NA":
        # parse a gromacs .top to a .itp file
        # print "ITP file created!"
        top.removeMolInfor(args.out)

if __name__ == "__main__" :
    # cd to PWD
    os.chdir(os.getcwd())

    if len(sys.argv) <= 2 :
        help_information = '''
        The tool is used for automating Gromacs based MD simulation.
        Several independent tools are provided.

        Separated tools including index, gentop, autoMD etc.

        Usage:

        python autoMD.py -h
        python autoMD.py gentop -h
        python autoMD.py index -h
        python autoMD.py extract -h
        '''
        print(help_information)
        sys.exit(1)
    else :
        if str(sys.argv[1]) == "index" :
            #print("Show help message: ")
            index = PdbIndex()
            index.genGMXIndex()
        #if 1 :
        elif str(sys.argv[1]) == "gentop" :
            runGenTop()

        elif str(sys.argv[1]) == "extract" :
            runExtract()
