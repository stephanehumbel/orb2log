#!/usr/bin/env python3
# put an orb file into a compatible gamess.log file for visualization in macmolplot
# usage  orb2log X.log Y.orb 
# generate YX.log that contains the orbs in proper format
# Note that the D orbitals are order correctly (difference between XMVB and GAMESS)
#
#
# take the MCSCF from 
## $MCSCF CISTEP=GUGA $END
## $DRT    GROUP=CS  FORS=.TRUE. NMCC=13 NDOC=2 NVAL=2 $END
## $GUGDIA NSTATE=2 $END
## $GUGDM2 WSTATE(1)=1.0,0.0 $END
## works also for RHF calculation (detects "EIGENVECTORS")
## dirty code to clean

import os
#import cclib
import sys
import re
#import routines
import numpy as np

def write_orb(filename, coeffs, indices, deb, fin):
#    print('write_orb',filename,len(coeffs),len(indices))
    if filename == 'screen':
        print('write_orb:',end=''  )
        for i in range(deb , fin):
            print(f"{len(indices[i]):4d}",end='')
            if (i+1) % 20 == 0:
                print() 
            # si on a plus de 50 valeurs, on saute une ligne
        print()
        for i in range(deb , fin):
            print(f"# ORBITAL {i+1:4d}  NAO = {len(indices[i]):4d}")
            count = 0
            for j in range(len(indices[i])):
                print(f"{float(coeffs[i][j]):13.10f}{(indices[i][j]):4d}  ",end='')
                #print(f"{float(coeffs[i][indices[i][j]]):13.10f}{(indices[i][j])+1:4d}  ",end='')
                count += 1
                if (j+1) % 4 == 0 and j != len(indices[i])-1:
                    print()
            print()
    else:        
        with open(filename, 'a') as f: 
            for i in range ( deb , fin):
                f.write(f"{len(indices[i]):4d}")
                if (i +1) % 20 == 0:
                         f.write("\n")

            f.write("\n")
            for i in range (deb , fin):
                f.write(f"# ORBITAL {i+1:4d}  NAO = {len(indices[i]):4d}      routines.write_orb({filename})\n")
                count = 0   
                for j in range(len(indices [i])):
#                    print(' coeffs(',i,',',j,')=',coeffs[i][1][indices[i][j]],end='')
                    #f.write(f"{coeffs[i][1][indices[i][j]]:13.10f}{(indices[i][j])+1:4d}  ")
                    f.write(f"{coeffs[i][j]:13.10f}{(indices[i][j]):6d}  ")
                    count += 1
                    if (j+1) % 4 == 0 and j != len(coeffs[i])-1:
                        f.write("\n")
                f.write("\n")
        f.close()

def print_lin_matrix(title,matrix):
    print('- ',title,end=' ')
    print('- mtrx printed in lines ------------')
    for i in range(len(matrix)):
        print(title,'[',f"{i:4d}",']',end="\t")
        for j in range(len(matrix[i])):
            print(f"{matrix[i][j]:7.3f}", end="\t")
            if (j + 1) % 10 == 0:
                print()
        print()

def print_matrix(title,matrix):
    print('- ',title,end=' ')
    print('- mtrx printed in lines ------------')
    for i in range(len(matrix)):
        print(title,'[',f"{i:4d}",']',end="\t")
        for j in range(len(matrix[i])):
            print(f"{matrix[i][j]:7.3f}", end="\t")
            if (j + 1) % 10 == 0:
                print()
        print()


def CS_weight(i,Ci,Sij):
    wi=Ci[i]*np.dot(Sij,Ci)[i]
#    print('CS_weight(',i,')=',wi)
    return wi
def Offset_conf(conf,offset):
    # offset the conf by the largest MO in VB conf
    # cares of the : in the confs
    # usefull to detect the max/min of the orbitals.
    k=0
    l=0
    inidi =[]
    new_indices=[]
    while True:
        if k >= len(conf):
            break
        else:
            l=0
        indices=re.split(' +|\n',conf[k])  
        #print(len(indices),'indices',indices,end=' ')
        llist=''
        while True:
            if l >= len(indices):
                break
            else:
                if indices[l] != '':
                    if ':' in  indices[l] :
                        indic=re.split(':',indices[l])  
                        debut=int(indic[0])+offset
                        fin=int(indic[1])+offset
                       # print('debut',debut,'fin',fin)
                        inidi.append(str(debut)+':' +str(fin)  )
                        llist+=str(debut)+':' +str(fin)+' '
                    else :
                        inidi.append(str(int(indices[l])+offset))
                        llist+=str(int(indices[l])+offset)+' '
            #print(inidi)
            l+=1
        new_indices.append(llist) 
        k+=1
    i=0
    tampon=[]
    if offset>=1:   
        str_ofset="1:"+str(offset)+'   '
        #print('str_ofset',str_ofset )
    else:
       str_ofset =" "
    while i < len(new_indices):
       str_ofset =" "
       tampon.append(str_ofset+new_indices[i])
       i+=1
    #print('new_indices',new_indices)
    return tampon

def Get_CIVECT(CAS_file_name, state):
    offset=4
    if state == -1:  #  xmo file
        offset=2   
        pos_CASvect,line=detect_keyword(CAS_file_name,"******  COEFFICIENTS OF STRUCTURES", 0)
    else:            # log file
        pos_CASvect,line=detect_keyword(CAS_file_name, "LAGRANGIAN CONVERGED", 0)
        search_state=""+str(state)+'  ENERGY'
        pos_CASvect,line=detect_keyword(CAS_file_name, search_state, pos_CASvect)
    if (pos_CASvect == -1):
        print("Get_CIVECT  no converged Ci coeffs in  ",CAS_file_name)
        print("for state :",state)
        quit(   )
    pos_CASvectfin = routines.detect_blank(CAS_file_name,pos_CASvect+offset)
    CI_SIZE=pos_CASvectfin-1-(pos_CASvect+offset-1)
    # on decale de 3 lignes pour .log et de 2 lignes pour xmo 
    CAS_conf, CAS_vect=Read_CIVECT(CAS_file_name, pos_CASvect+offset-1,pos_CASvectfin-1,offset)
    return CAS_conf,CAS_vect

def est_reel(chaine):
    return '.' in chaine

 # check the lenght & size
def check_size(file_name,pos,fin):
    with open(file_name, 'r') as file:
        k=0 # line counter
        noa=0
        nom=0
        NOA=-1  
        NOM=0
        #print('--')
        while k < pos: # skip the first pos lines
            file.readline() 
            k+=1
        newblock=False
        while k < fin: # will stop 
            line = file.readline() 
            if len(line.split()) >=1:
                noa+=1
                newblock=True
                nom=len(line.split())-1
                #print('',len(line.split()),end=' ')
            elif newblock and len(line.split()) == 0:
                if NOA==-1:
                    NOA=noa-1
                    noa=0
#                    print('NOA=',NOA,end=' ')
                NOM=NOM+nom
#                print('size=',NOA,'x',NOM,end=' ')
                newblock=False
            k+=1
#        print('NOM=',NOM,end=' ')
    file.close()
    return NOA,NOM

def read5cols(file_name,length,size,pos,fin):# read a file with 3 blank lines, +1 to skip, then blocks of columns of length lines
                                         #  reading starts at pos. end reading when size vectors are read    
                                         # pos  A                 ******  OVERLAP OF VB STRUCTURES  ******
                                         # .    B
                                         #      C
                                         #      D           1            2            3            4            5
                                         #       1       1.000000     0.283693     0.222389    -0.370804    -0.341846
                                         #       2       0.283693     1.000000     0.283693    -0.149626    -0.132904
                                         #       3       0.222389     0.283693     1.000000    -0.370806    -0.122930
                                         #       4      -0.370804    -0.149626    -0.370806     1.000000     0.262162
    if size==-1:
        # should be able to read undetermined size and length
        # and the end of the section 
        size=length
        #continue
  #  nblock= (fin-pos -5)/length
  #  print('nblock= environ',nblock, (fin-pos -5-3*(nblock-1))/length, (fin-pos -5-3*(nblock-1))//length)
  #  print('size= environ',  (fin-pos -5-3*(nblock-1))/nblock)
    Stot=np.zeros((length,size))
    ioa=0
    nblock=0
    line='blanck'
    space=0
    firstom=0
    lenlue=0
    with open(file_name, 'r') as file:
        k=0 # line counter
        print('--')
        while k < pos: # skip the first pos lines
            file.readline() 
            k+=1
        while k < fin: # will stop 
            line = file.readline() 
            k+=1
            if len(line.split()) == 0: #line is empty=space between block
                if space==0: 
    #                print('new block starts, nblock',nblock,end=' ')
                    firstom=firstom+lenlue
                    ioa=0
                line = file.readline() 
                k+=1
                space+=1
            else: 
                if space !=0: # new block
                    space=0
                    nblock+=1
                    line = file.readline() 
                    k+=1
                lenlue=len(line.split() )-1
                for i, val in enumerate(line.split()[1:]):
                    iom=i+firstom
    #                print('S(',iom,',',ioa,')','=',val,end=' ')
                    Stot[iom][ioa] = float(val)
                ioa+=1
        print("| ------",iom+1,' columns have been read')
        file.close()
        return Stot


def read5OM_LOG(file_name,nOM,nOA,pos,fin):# read a file with 3 blank lines, +1 to skip, then blocks of columns of length lines
       
#          ------------------------
#          MCSCF OPTIMIZED ORBITALS
#POS          ------------------------
#A   
#B                      1          2          3          4          5
#C                  -20.5795   -11.3211   -11.2497   -11.2447    -1.3983
#D                     A'         A'         A'         A'         A'
#->  1  C  1  S    0.000006  -0.000049   0.982745  -0.159416  -0.003006
#    2  C  1  S    0.000037  -0.000252   0.027229  -0.004718   0.007837
#    3  C  1  X   -0.000020   0.000012  -0.000066   0.000220  -0.003376
#    4  C  1  Y   -0.000016   0.000146  -0.000088   0.000133  -0.004117
#    5  C  1  Z   -0.000000   0.000000   0.000000   0.000000   0.000000
#    6  C  1  S   -0.000085  -0.000083  -0.010422   0.006291   0.005935
#    7  C  1  X    0.000035  -0.000036   0.000353  -0.001412  -0.002224

    if nOM==-1:
        # should be able to read undetermined size and length
        # and the end of the section 
        nOM=nOA
        if nOM<5:
            nOM=5

        #continue
  #  nblock= (fin-pos -5)/length
  #  print('nblock= environ',nblock, (fin-pos -5-3*(nblock-1))/length, (fin-pos -5-3*(nblock-1))//length)
  #  print('size= environ',  (fin-pos -5-3*(nblock-1))/nblock)

    #print('MO_CAS',nOM,nOA, pos,fin)
    #truc sale pour lire moins que 5 orbitales, 
    if nOM<5:
            nOM=5

  
    MO_CAS=np.zeros((nOM+1,nOA+1))
    ioa=0
    iom=0
    nblock=0
    line='blanck'
    space=0
    firstom=0
    lenlue=0
    left_to_read=nOM
    next_block_size=5
    readentete=True
    entetes=[]
    with open(file_name, 'r') as file:
        k=0 # line counter
        #print('--')
        while k < pos: # skip the first pos lines
            line=file.readline() 
#            print('MO_CAS',line[1:20], k,fin)
            k+=1
        #print('MO_CAS',line[1:20], k,fin)
        while k < fin or left_to_read >0 : # will stop 
            line = file.readline()
            #print('MO_CAS',line, k)
            k+=1
            if len(line.split()) == 5: #line is the header of a block
                continue
            if len(line.split()) == 0: #line is empty=space between block
                if space==0:                        
                    line = file.readline() 
                    #print('MO new block starts, nblock',nblock,end=' ')
                    firstom=firstom+lenlue
                    left_to_read=nOM-firstom
                    next_block_size=left_to_read % 5
                    if left_to_read > 5: # blocks are normaly 5 colums
                        next_block_size=5
                    if left_to_read < 1:
                        break
                    #print('MO new block starts, nblock',nblock,left_to_read,next_block_size,end=' ')
                    ioa=0
                #print('-')
                k+=1
                space+=1
            else: 
                if readentete:
                    entete=line[0:14]
                    entetes.append(entete)
                  #  print('entete',nOA,entete,'|')
                    if len(entetes)>=nOA:
                        readentete=False
                    continue

                if space !=0: # new block
                    space=0
                    nblock+=1
                   # line = file.readline() 
                    k+=1
                lenlue=len(line.split() )-4
                for i, val in enumerate(line.split()[4:4+next_block_size]):
                    iom=i+firstom
                   # if left_to_read < 5:
                   #     print('RRR',line.split()[4:],end=' ')
                   # print('|%|',nOM,f"{float(val):7.3f}",end=' ')
                    #print(f"{float(val):7.3f}))",end=' ')
                    #print('MO(',iom,',',ioa,')','=',val,end=' ')
                    MO_CAS[iom][ioa] = float(val)
                ioa+=1
                #print(ioa,end=' ')
        print("| ------",iom+1,' columns have been read')
        file.close()
        return entetes,MO_CAS


def Read_CIVECT(file, pos,fin,offset):
    #offset=4 pour .log et 2 pour .xmo
    # indique la position a lire
    with open(file, 'r') as f:
        lines = f.readlines()[pos:fin]
        n=fin-pos
#        print('n',pos,n,fin)
        CAS_vect = np.zeros(n)
        CAS_conf = []
        for i in range(n):
# Det,       CAS_vect[i] = float(re.split('\|+| +|\n', lines[i])[7])
# GUGA
            CAS_vect[i] = float(re.split(' +',lines[i])[2])
            if offset == 4: # log
                str=re.split(' +|\n',lines[i])[3]
                CAS_conf.append(str)
            else: # xmo
                #CAS_conf[i] = re.split(' +|\n',lines[i])[3:]
                str=''
                k=4
                while True:
                    if re.split(' +|\n',lines[i])[k]!= '':
                        str+=' '+re.split(' +|\n',lines[i])[k]   
#                        print(str,end=' ')
                    else:
                        break
                    k+=1
                CAS_conf.append(str)
                continue
    f.close()
    return CAS_conf, CAS_vect
def ls_dir(beginning,end,k):
    files = [file for file in os.listdir() if file.endswith(end) and file .startswith(beginning)]
    for file in files:
        print(file, end=' ')
        k+=1
        if k%5 == 0:
            print() 
    return  k

def reverse_reorder_OA(tab): # 
    reord = list(range(len(tab)))  # initial indices
    for i in range(len(tab)):
        if tab[i] == "XX" :
              #reord the 6 AO's 'XX', 'XY', 'XZ', 'YY', 'YZ', 'ZZ'
              # into            'XX', 'YY', 'ZZ', 'XY', 'XZ', 'YZ'
                              #   0     1     2     3     4     5
                              #   0     3     5     1     2     4
            reord[i + 1:i + 6] = [i + 3, i + 5, i + 1, i + 2, i + 4]
            i+=5
    return reord       
               
def Read_INT(line, keyword):
    words = re.split(' +|=+|\n',line)
    strin=' '.join(words) # twice to remove all the empty strings 
    words = re.split(' +|=+|\n',strin)
    index = words.index(keyword)
    inte = int(words[index + 1])
    return inte

def detect_keyword(file_path, keyword, start_line):
    '''
    Description : get the line number of the *last* occurence of keyword in file_path.
    Args:
        parameter (file_path): name of the file (possibly the path?).
        parameter (keyword): string to find
        parameter (start_line): the line the search starts from

    Returns:
        type: returns the line number (lin_num) of the *last* occurence of keyword or -1 if keyword is not found
    '''
#    print('detect_keyword',file_path,keyword,start_line)
    with open(file_path, 'r') as file:
        ret_num = -1                # not found
        ret_line=''                      # not found
        for line_num, line in enumerate(file, 1):  
            if line_num >= start_line and keyword in line:
                ret_line = line
                ret_num = line_num
    return ret_num, ret_line 


def detect_next_keyword(file_path, keyword, start_line):
    '''
    Description : get the line number of the *next* occurence of keyword in file_path.
    Args:
        parameter (file_path): name of the file (possibly the path?).
        parameter (keyword): string to find
        parameter (start_line): the line the search starts from

    Returns:
        type: returns the line number (lin_num) of the *next* occurence of keyword or -1 if keyword is not found
    '''
#    print('detect_keyword',file_path,keyword,start_line)
    with open(file_path, 'r') as file:
        ret_num = -1                # not found
        ret_line=''                      # not found
        for line_num, line in enumerate(file, 1):  
            if line_num >= start_line and keyword in line:
                ret_line = line
                ret_num = line_num
                break
    return ret_num, ret_line 



def is_pi(orb):
    num_zeros = orb.count(0)
    ratio = num_zeros / len(orb)
    return ratio > 0.25

def read_orb(file_name):
    '''
    Description :  Reads coef and ao in  an .orb file 
    Args:
        (file_name):  the file .orb

    Returns:
        the tables all_coeffs and all_aos
    '''
    all_coeffs = []  
    all_aos = []     
    coef = []  
    ao = []
    om1 = False  
    with open(file_name, 'r') as file:
        for line in file:
            line = line.strip()
            if ("#" in line):
                if om1:  
                    all_coeffs.append(coef)
                    all_aos.append(ao)
                    coef = []  
                    ao = []
                else :
                    om1 = True
            elif om1 :
                values = line.split()
                for i in range(0, len(values), 2):
                    coef.append(float(values[i])) # add the last to vectors
                    ao.append(int(values[i + 1]))

        # last orb must be updated        
        all_coeffs.append(coef)
        all_aos.append(ao)
    return all_coeffs, all_aos
def reorder_OA(tab):
    reord = list(range(len(tab)))  # initial indices
    for i in range(len(tab)):
        if tab[i] == "XX" :
              #reord the 6 AO's 'XX', 'YY', 'ZZ', 'XY', 'XZ', 'YZ'
              # into            'XX', 'XY', 'XZ', 'YY', 'YZ', 'ZZ'
                               #  0     1     2     3     4     5
                               #  0     3     4     1     5     2
            reord[i + 1:i + 6] = [i + 3, i + 4, i + 1, i + 5, i + 2]
            i+=5
    return reord       

def make_vec(coeffs,indices,nao):
    ''' from short .orb, Returns  long vectors of MOs with all the zero AOs '''
    #print("make_vec",len(coeffs),coeffs,indices,nao)
    long_ao = []
    long_coeffs = []
    ao = list(range(1, nao + 1))
    coeffsis = [0.0] * nao
    #print("@@@@@@@@",len(coeffs),nao,indices)
    #print("------")
    for i in range(0,len(coeffs)):
        coeffsis = [0.0] * nao
        #print("i=",i," OM n°",i+1, coeffs[i], end='\n ')
        for j in range(0,len(coeffs[i])):
            #print("j",j, coeffs[i][j],indices[i][j], end='')
            coeffsis[indices[i][j]-1] = coeffs[i][j]
            #print(coeffsis,"\n")
        long_ao.append(ao)
        long_coeffs.append(coeffsis)
    #print("FFFff", long_coeffs,"\n et fin")
    return long_coeffs,long_ao

def write_orbs_as_gamess(orb, entete, filein, fileout, debut, fin):
    # Ouvrir le fichier d'entrée et lire les lignes
    with open(filein, 'r') as f:
        lines = f.readlines()

    # Extraire l'en-tête et le pied de page
    header = lines[:debut]
    footer = lines[fin:]

    # Ouvrir le fichier de sortie en mode écriture
    with open(fileout, 'w') as file:
        # Écrire l'en-tête
        file.writelines(header)
        norb=len(orb)
        # Parcourir les orbitales par blocs de 5 colonnes
        i=0 
        while  i <=(len(orb)):
            a_ecrire=max(min(5,norb-i),(len(orb)-i)%5)
            file.write(f"             ") 
            for j in range(a_ecrire):
                file.write(f"{i+j+1:10} ")        
            file.write("\n")        
            
            file.write(f"               ") 
            for j in range(a_ecrire):
                file.write(f"{float(j+i+1):10.4} ")        
            file.write("\n")     
            file.write(f"                     ") 
            # Écrire les labels des orbitales 
            pi_labels=[]
            for j in range(a_ecrire):
                if is_pi(orb[j+i]):
                    pi_label="Loc"
                else:
                    pi_label="Del"
                pi_labels.append(pi_label)
            #    print(pi_labels,end=' ')
            for j in range(a_ecrire):
                file.write(f"{pi_labels[j]}{j+i+1:<8}")       
            file.write("\n")           
          
 #           print('a_ecrire',a_ecrire)
            # Écrire l'en-tête de l'orbitale
            for iao in range(len(orb[i])):
                file.write(f"{entete[iao]}  ")
            # Écrire les valeurs des orbitales pour cette ligne
                for j in range(a_ecrire):
                    if j < len(orb[i]):
                        file.write(f"{orb[j+i][iao]:10.6f} ")
                    else:
                        file.write("           ")
                file.write("\n")        

  #         for j in range(len(orb)):
  #              file.write(f"{orb[j][i]:10.6f} ")
            
            # Sauter à la ligne après chaque orbitale
            file.write("\n")
            
            # Sauter 3 lignes après chaque bloc de 5 colonnes
            if (i + 1) % 5 == 0:
                file.write("#@#\n\n\n")
            i+=5
        # Écrire le pied de page
        file.writelines(footer)

# read types in basis
def read_basis(file_name):
    num, _ = detect_keyword(file_name, "OPTIMIZED ORBITALS", 0)
    types = []
    if num == -1:
         num, _ = detect_keyword(file_name, "EIGENVECTORS", 0)
    print('ReadBasisnum=',num)
    with open(file_name, 'r') as file:
        for line_num, line in enumerate(file, start=1):
            if line_num >= num + 6:
                if line.strip() == '':
                    break
                values = re.split(r'\s+', line.strip())
                types.append(values[3])

    return types
# =========------------------------------------------------------
# == main =------------------------------------------------------
# =========------------------------------------------------------
print("+--------copyorb.py - SH 2025 ---------------------")
print('| copyorb.py file1CAS.log file2.orb ',len(sys.argv), 'arguments')
print("+- must have .log and a .orb              ---")
print("+--------                     ---------------------")
vect=[]
entete=[]
if len(sys.argv) == 1:
    print("+- must have .log and a .orb !!              ---")
    quit()
if len(sys.argv) == 2:
#For an xmo file')
    CAS_file = sys.argv[1]
    CAS_file_name, CAS_file_ext = os.path.splitext(CAS_file)
#For a  GAMESS log file')
    if os.path.exists(CAS_file_name+'.log'):
        print (CAS_file_name+'.log exists must be with these types of keywords \n $MCSCF CISTEP=GUGA MAXIT=200 QUAD=.F. $END\n $GUGDIA NSTATE=2 $END $GUGDM2 WSTATE(1)=1.0,0.0 $END')
        log_file=CAS_file_name+'.log'
        if os.path.exists(CAS_file_name+'.orb'):
            orbfile=CAS_file_name+'.orb'
            print('cool, ',orbfile,'  exists ')
        else:    
            print('please provide a .orb file ')
            orbfile=input('give the .orb file:')
            if not(os.path.exists(orbfile)):
                print('### >>> ',orbfile,'  not found <<<<< ##')
                quit()  
    else:
        quit()
if len(sys.argv) == 3:
    log_file = sys.argv[1]
    orbfile=sys.argv[2]

    if True: # pour eviter de de-indenter 
        #CAS_conf,CAS_vect=Get_CIVECT(log_file,1)
        #toprint=routines.make_conf_from_gamess(CAS_conf)
        #num,line=detect_keyword(log_file, "SPIN MULTIPLICITY", 0)    
        #MULT=Read_INT(line,"MULTIPLICITY")
        num,line=detect_keyword(log_file, "NUMBER OF CARTESIAN GAUSSIAN", 0)  
        NBASIS=Read_INT(line,"FUNCTIONS")
        num,line=detect_keyword(log_file, "NUMBER OF OCCUPIED ORBITALS (ALPHA)", 0)
        norb=Read_INT(line,"(ALPHA)")
        type_OA=read_basis(log_file)
        print(norb, len(type_OA))
        
        reord_OA=[]
        #reord_OA=routines.reorder_OA(type_OA)
        print(type_OA,reord_OA)
        #new_vect=[]
        pos,line=detect_keyword(log_file, "MCSCF OPTIMIZED ORBITALS", 0)
        if pos == -1:
            pos,line=detect_keyword(log_file, "MOLECULAR ORBITALS", 0)
        print("pos=",pos)
        posfin,line=detect_next_keyword(log_file,"END OF",pos)
        if posfin == -1:
            posfin=line=detect_next_keyword(log_file,"CALCULATION",pos)
        pos+=2
        print('  ;;; read files : ',norb,log_file,pos,posfin,end=':')
           
        entete,vect=read5OM_LOG(log_file,norb,NBASIS,pos,posfin)
       # reord_OAv=reorder_OA(type_OA)
       ## new_vect=[]
       # indices=[]
       # for j in range(len(vect)):
       #  new_orb=[]
       ##  indic=[]
       #  for i in range(len(reord_OA)):
       #         indic.append(i)
       #         new_orb.append(vect[j][reord_OA[i]])
       #  new_vect.append(new_orb)
       #  indices.append(indic)
       # print(vect)
       # write_orb("screen", new_vect, indices, 0, len(vect))
 #
    print('')
    new_vvect=[]
# read the orbfile  
    orb_coeffs,orb_aos=read_orb(orbfile)
    long_coeff,long_ao=make_vec(orb_coeffs,orb_aos,len(type_OA))
# reordonne a cause des orbitales d de XMVB qui ne sont pas dans le meme ordre que dans GAMESS
    reord_OA=reverse_reorder_OA(type_OA)
    print(reord_OA)
    for j in range(len(long_coeff)):
        new_orb=[]
        for i in range(len(reord_OA)):
                new_orb.append(long_coeff[j][reord_OA[i]])
        new_vvect.append(new_orb)

# choisi l'output

    log_n, log_e = os.path.splitext(log_file)
    orb_n, orb_e = os.path.splitext(orbfile)
    logout=orb_n+log_n+'.log'
    print('|  I wrote the new orbfile in ',logout)
    write_orbs_as_gamess(new_vvect,entete,log_file,logout,pos,posfin)


    quit()
quit()
#------------------------------------------------------------------------#
#------------------------------------------------------------------------#

