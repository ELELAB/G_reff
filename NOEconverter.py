#!/usr/bin/env python 
#Dependencies
import sys
import itertools


#Classes
#An automated nested dictionary class makes
class Vividict(dict):
    def __missing__(self, key):
        value = self[key] = type(self)()
        return value

#Functions
#Reads information of the protein from pdb file, sequence, residues and atoms outputs in a nested dict as dict[residue][protein]=[list of atoms] 
def readprotein(pdbfile,convfile,shift=0,sequence=Vividict(),rawsequence=Vividict()):
    f = open(pdbfile,"r")
    converter = readconversiontable(convfile)
    for line in f:
        split = line.split()
        if split[0]=="ATOM":
            atom = converter[split[3]][split[2]]
            #print atom,split[3],split[2] 
            #print atom, split[2],int(split[5])+shift
            sequence[int(split[5])+shift].setdefault(split[3],[]).append(atom) #Number,Restype,atoms
            rawsequence[int(split[5])+shift].setdefault(split[3],[]).append(split[2]) #Number,Restype,atoms
        elif split[0]=="ENDMDL":
           break
    f.close()
    return sequence, rawsequence

#Reads the conversion table file in order to convert from atom names given in pdb file to NOE file
def readconversiontable(tablefile):
    f = open(tablefile,"r")
    convtable = Vividict()
    for line in f:
        split = line.split()
        convtable[split[0]][split[2]]=split[1]    
    f.close()
    return convtable

#Reads NOE file and outputs the converted file 
def openNOE(distances,protein,rawprotein):
    f = open(distances,"r")
    f2 = open("NOE_converted.dat","w")
#    for idx,line in enumerate(f):
#        split=line.split()
#        if int(split[2]) > 300:
#            resid1=int(split[2])-266
#            chain1="B"
#        else:
#            resid1=int(split[2])+34
#            chain1="A"
#        if int(split[8]) > 300:
#            resid2=int(split[8])-266
#            chain2="B"
#        else:
#            resid2=int(split[8])+34
#            chain2="A"
#        atomtype1=testatom(split[5],resid1,protein,rawprotein)
#        atomtype2=testatom(split[11],resid2,protein,rawprotein)
#        print atomtype1"""
    for idx,line in enumerate(f):
        split =line.split()
        if split[0]=="assign":
            #Renaming split for convinience sake
            split = line.split()
            #Renumbering of residues and assignment of chains
            if int(split[2]) > 300:
                resid1=int(split[2])-266
                chain1="B"
            else:
                resid1=int(split[2])+34
                chain1="A"
            if int(split[8]) > 300:
                resid2=int(split[8])-266
                chain2="B"
            else:
                resid2=int(split[8])+34
                chain2="A"
            atomtype1=testatom(split[5],resid1,protein,rawprotein)
            atomtype2=testatom(split[11],resid2,protein,rawprotein)
            classtype=idx #for now
            reff_exp=split[15]
            upper=0
            lower=0
            firstline=True
            for prod in itertools.product(atomtype1,atomtype2): #Itertools makes all possible combination of list1 and list2
                if firstline==True: #First line is different
                    f2.write("%s %s %s %s %s %s %s %s %s %s\n"%(resid1,prod[0],chain1,resid2,prod[1],chain2,classtype,reff_exp,upper,lower))
                    firstline=False
                elif firstline==False:
                    f2.write("%s %s %s %s %s %s %s\n"%(resid1,prod[0],chain1,resid2,prod[1],chain2,classtype))
    f.close()
    f2.close()

#Taking atom in protein and converting to list, if # is in BMRB - all cooresponding atoms are found
def testatom(atomstring,resID,protein,rawprotein):
    if resID not in protein:
        print resID, "residue missing in pdb"
        return [atomstring]
    elif "#" in atomstring: #Check if # is instead of one number
        atomlist=[]
        for n in range(1,4):
            atomreplace=atomstring.replace("#",str(n))
            for idx,atom in enumerate(protein[resID][next(iter(protein[resID].keys()))]):
                if atomreplace == atom:
                    if resID == 104:
                        print "s" +atomstring+str(resID)
                        print atom
                        print protein[resID].keys()
                    atomlist.append(rawprotein[resID][next(iter(rawprotein[resID]))][idx])
        if len(atomlist)<1: #Check if # is instead of two numbers
            for n in range(1,4):
                for j in range(1,4):
                    atomdreplace=atomstring.replace("#",str(n)+str(j))
                    for idx,atom in enumerate(protein[resID][next(iter(protein[resID].keys()))]):
                        if atomdreplace == atom:
                            atomlist.append(rawprotein[resID][next(iter(rawprotein[resID]))][idx])
        return atomlist
    else:
        for idx,atom in enumerate(protein[resID][next(iter(protein[resID].keys()))]):
            if atomstring == atom:
                atomstring2=rawprotein[resID][next(iter(rawprotein[resID]))][idx]
                return [atomstring2]
        print "Error", atomstring, resID

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Converts BMRB NOE to g_reff format')
    parser.add_argument("-pdb",help="pdb structure",required=True)
    parser.add_argument("-noe",help="NOE_BMRB file",required=True)
    parser.add_argument("-c", help="Atom conversion table",required=True)

    args = parser.parse_args()
    protein,rawprotein=readprotein(pdbfile=args.pdb,convfile=args.c)
    openNOE(args.noe,protein,rawprotein)
    
