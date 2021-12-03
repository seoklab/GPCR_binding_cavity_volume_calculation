#!/usr/bin/env python
import numpy as np
import sys
import torch
import os,sys
import copy
from time import time
import argparse
#########################################################################
# parameter set up
torch.set_num_threads(1)
radius_dict = {'H':1.2, 'C':1.77, 'N':1.66, 'O':1.50, 'F':1.46, 'Cl':1.82, 'Br':1.86, 'I':2.40, 'P':1.90, 'S':1.89}
element_list = ['C','N','O','S']
n_elem=81; grid_width = 0.5; default=(0-((n_elem-1)*grid_width)/2.0)
##########################################################################
class Atom:
    def __init__(self):
        pass
    def parse_atom_line(self, line):
        self.header = line[:6]
        self.atm_no = int(line[6:11])
        self.atm_name = line[12:16]
        self.red=line[16]
        self.res_name = line[17:20].strip()
        self.chain_id = line[21]
        self.res_no = int(line[22:26])
        self.R = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
class PDB:
    def __init__(self, fn):
        self.atom_s = []
        self.read_pdb(fn)
    def read_pdb(self, fn):
        with open(fn) as fp:
            for line in fp:
                if not line.startswith('ATOM'):
                    continue
                atom = Atom()
                atom.parse_atom_line(line)
                self.atom_s.append(atom)
    def write_pdb(self, out_fn):
        wrt=[] 
        for atom in self.atom_s:
            line = '%-6s%5i %-4s %3s %s%4i    %8.3f%8.3f%8.3f\n'\
                    %(atom.header, atom.atm_no, atom.atm_name, atom.res_name, atom.chain_id,\
                    atom.res_no, atom.R[0], atom.R[1], atom.R[2] )
            wrt.append(line)
        with open(out_fn,'wt')as fp:
            fp.writelines(wrt)

def trs_to_cntr(prot,trp_res_no,cut_Nterm,cut_Cterm):
    # Translate protein to set center of protein 0,0,0
    new_atom_s=[];num_atom=0;val=np.array([0.0,0.0,0.0])
    for atom in prot.atom_s:
        if (atom.res_no <cut_Nterm) or (atom.res_no>cut_Cterm):
            continue
        num_atom+=1;val+=atom.R
    val= val/float(num_atom)
    for atom in prot.atom_s:
        atom.R=atom.R-val
        new_atom_s.append(atom)
    prot.atom_s=new_atom_s # update coordinates of protein
    # Translate protein to set z-coordinate of C-alpha atom of toggle switch TRP 
    new_atom_s=[]
    for atom in prot.atom_s:
        if atom.res_no ==trp_res_no and atom.atm_name.strip()=='CA':
            ref=copy.deepcopy(atom.R)
    for atom in prot.atom_s:
        atom.R[-1]=atom.R[-1]-ref[-1]
        new_atom_s.append(atom)
    prot.atom_s=new_atom_s
    #
    return prot

def truncate_side(mat,cnt_s,exclude_tm1):
    if exclude_tm1:
        print ("TM1 is excluded for side truncation")
        cnt_s=cnt_s[1:7]
    for k in range(n_elem): # iteration on z coordinate
        wrt=[]
        corner_s=call_corner_s(cnt_s,-k)
        x_range_s=[]
        for x in range(len(corner_s)):
            if x==len(corner_s)-1:
                x_range_s.append(sorted([corner_s[x][0],corner_s[0][0]]))
                continue
            x_range_s.append(sorted([corner_s[x][0],corner_s[x+1][0]]))
        line_s=[]
        for x in range(len(corner_s)):
            if x==len(corner_s)-1:
                line_s.append(eqn_line(corner_s[x],corner_s[0]))
                continue
            line_s.append(eqn_line(corner_s[x],corner_s[x+1]))
        for i in range(n_elem):
            for j in range(n_elem):
                grid_crd=g_crd(i,j,k)
                r_stat=identify_inout(line_s,grid_crd,x_range_s)
                if not r_stat :
                    mat[i,j,k,0]=2
    return mat



    
def truncate_roof(mat,tip_s):
    # Define 5 planes with tip points of following  1,2,7 / 2,7,3 / 3,7,4 / 7,4,6 / 4,6,5. Arbitraily defined....
    pre_line_s=[[2,7],[7,3],[7,4],[4,6]]
    line_s=[]
    for x in pre_line_s:
        line_s.append(eqn_line(tip_s[x[0]-1],tip_s[x[1]-1]))
    pre_plane_s=[[1,2,7],[2,3,7],[3,4,7],[4,6,7],[4,5,6]]
    plane_s=[]
    for x in pre_plane_s:
        plane_s.append(eqn_plane(tip_s[x[0]-1],tip_s[x[1]-1],tip_s[x[2]-1]))
    wrt=[]
    for i in range(n_elem):
        for j in range(n_elem):
            for k in range(n_elem):
                if mat[i,j,k,0]==2:
                    continue
                grid_crd=g_crd(i,j,k)
                plane_num=identify_plane_num(grid_crd,line_s)
                plane=plane_s[plane_num]
                if calc_plane(grid_crd,plane) >0:
                    mat[i,j,k,0]=2
                else:
                    wrt.append('ATOM      0  H1  DUM L        %8.3f%8.3f%8.3f\n'%(grid_crd[0],grid_crd[1],grid_crd[2]))
    return mat,wrt

def eqn_line(p1,p2):
    x1,y1,x2,y2=p1[0],p1[1],p2[0],p2[1]
    slope=(y2-y1)/(x2-x1)
    y_intercept=y1-slope*x1
    if slope>0:
        sign=1
    else:
        sign=-1
    return [slope,y_intercept,sign]
def eqn_plane(p1,p2,p3):
    #get plane for given 3 points
    v1=p2-p1;  v2=p3-p1;
    vn=np.cross(v1,v2)
    const=-1*np.dot(vn,p1)
    return vn,const
##TODO
def call_corner_s(cnt_s,z_coord):
    corner_s=[[-9999,-9999,-9999] for x in range(len(cnt_s))]
    for x in range(len(cnt_s)):
        tm=cnt_s[x]; memo=[]
        for cnt in range(len(tm)):
            memo.append([np.linalg.norm(tm[cnt][-1]-z_coord*grid_width),cnt])
        memo.sort()
        corner_s[x]=tm[memo[0][1]]
    return corner_s

def identify_inout(line_s,grid_crd,x_range_s): 
    # check whether the given grid point is belong to inside of GPCR or outside of GPCR
    cand_s=[]
    for y,x_range in enumerate(x_range_s): 
        if x_range[0]<= grid_crd[0] <=x_range[1]:
            cand_s.append(y)
    n_contact=0
    memo=[]
    for x in cand_s:
        tmp=line_s[x]
        cont_y=tmp[0]*grid_crd[0]+tmp[1]
        memo.append(cont_y)
        if cont_y >grid_crd[1]:
            n_contact+=1
    if n_contact%2 == 0 :
        return False
    else:
        return True
def get_tips(cnt_s):
    return [ x[-1] for x in cnt_s ]
def calc_helix_center(prot,point_s):
    helix_resno_s=[]
    #define residues for calculating center of transmembrane helix
    for tm_idx,tip_resno in enumerate(point_s):
        if (-1)**tm_idx >0:
            helix_resno_s.append(range(point_s[tm_idx],point_s[tm_idx]+20))
        else:
            helix_resno_s.append(range(point_s[tm_idx]-19,point_s[tm_idx]+1))
    helix_coord_s =[[] for x in range(len(point_s))]
    for atom in prot.atom_s:
        if not atom.atm_name.strip()=='CA': # use only C-alpha coordinate
            continue
        for tm_idx,resno_range in enumerate(helix_resno_s):
            if atom.res_no in resno_range:
                helix_coord_s[tm_idx].append(atom.R)
    for x,_ in enumerate(helix_coord_s):
        helix_coord_s[x]=sorted(helix_coord_s[x] , key=lambda y: y[-1],reverse=True)
    cnt_s=[ [] for x in range(len(point_s))] # traces of helix center
    for tm_idx, _ in enumerate(helix_coord_s):
        for j in range(len(helix_coord_s[tm_idx])-3):
            cnt_s[tm_idx].append(np.mean(helix_coord_s[tm_idx][j:j+4],axis=0))
        for x in range(3,0,-1):
            cnt_s[tm_idx].append(np.mean(helix_coord_s[tm_idx][-x:],axis=0))
    # interpolation twice
    for _ in range(2):
        new_cnt_s=[]
        for x in cnt_s:
            tmp=[]
            for y in range(len(x)-1):
                tmp.append(x[y])
                tmp.append(np.mean([x[y],x[y+1]],axis=0))
            tmp.append(x[-1])
            new_cnt_s.append(tmp)
        cnt_s=new_cnt_s
    #sort
    for tm_idx, _ in enumerate(cnt_s):
        cnt_s[tm_idx]=sorted(cnt_s[tm_idx] , key=lambda y: y[-1],reverse=True)
    # center of tip point
    grid_crd=np.mean([new_cnt_s[x][-1] for x in range(7)],axis=0)
    wrt=[]; 
    wrt.append('ATOM      0  O1  DUM L        %8.3f%8.3f%8.3f\n'%(grid_crd[0],grid_crd[1],grid_crd[2]))
    # record helix center trace
    for x in cnt_s:
        for y in x:
            grid_crd=y
            wrt.append('ATOM      0  H1  DUM L        %8.3f%8.3f%8.3f\n'%(grid_crd[0],grid_crd[1],grid_crd[2]))
    return cnt_s,wrt

def calc_line(grid_crd,tmp):
    slope=tmp[0]
    work=tmp[1]
    return tmp[2]*(grid_crd[1] - (slope*grid_crd[0] + work))
def calc_plane(grid_crd,tmp):
    val=0.0
    for x in range(3):
        val+=tmp[0][x]*grid_crd[x]
    val+=tmp[-1]
    return val
def identify_plane_num(grid_crd,line_s):
    if calc_line(grid_crd,line_s[0]) >0:
        return 0
    elif calc_line(grid_crd,line_s[0])<0 and calc_line(grid_crd,line_s[1])>0:
        return 1
    elif calc_line(grid_crd,line_s[1])<0 and calc_line(grid_crd,line_s[2])>0:
        return 2
    elif calc_line(grid_crd,line_s[2])<0 and calc_line(grid_crd,line_s[3])>0:
        return 3
    elif calc_line(grid_crd,line_s[3])<0:
        return 4
def g_crd(i,j,k):
    return torch.FloatTensor([default+grid_width*i,default+ grid_width*j,-grid_width*k])
def gen_mat():
    return torch.zeros(n_elem,n_elem,n_elem,1)
def check_contact(res_R, cutoff_list, crd):
    contact = False; occupied = False
    crd_expand = crd.expand(res_R.size()[0], -1)
    distances = res_R - crd_expand
    distances = distances*distances
    distances = torch.sum(distances, 1)
    distances = torch.sqrt(distances)
    crit = distances - cutoff_list
    min_crit, index = torch.min(crit,0)
    if (min_crit > -grid_width/2.0 and min_crit < grid_width/2.0):
        contact = True
    if min_crit < 0:
        occupied = True
    return contact, occupied

def get_element(atom_name):
    if not atom_name[0].isdigit():
        element = atom_name[0]
    else:
        element = atom_name[1]
    return element

def check_volume_and_sasa(bsite_residues, cntr_crd,mat):
    res_R, cutoff_list = get_R_in_bsite(bsite_residues)
    min_crd = torch.FloatTensor([x-((n_elem-1)*grid_width)/2.0 for x in cntr_crd])
    n_unoccupied = 0; wrt=[]
    for i in range(n_elem):
        for j in range(n_elem):
            for k in range(n_elem):
                if mat[i][j][k][0]==2:
                    continue
                grid_crd = torch.FloatTensor([default+grid_width*i, default+grid_width*j, -grid_width*k])
                contact, occupied= check_contact(res_R, cutoff_list, grid_crd)
                if not occupied:
                    n_unoccupied += 1
                    wrt.append('ATOM      0  H1  DUM L        %8.3f%8.3f%8.3f\n'%(grid_crd[0],grid_crd[1],grid_crd[2]))
    dat_wrt=[]
    dat_wrt.append("Number of unoccupied grid : %i\n"%n_unoccupied)
    dat_wrt.append("Single grid point represent %5.3f A^3\n"%(grid_width **3))
    dat_wrt.append("Volume : %5.3f A^3\n"%(n_unoccupied*(grid_width**3)))
    print ("Number of unoccupied grid : %i"%n_unoccupied)
    print ("Single grid point represent %5.3f A^3"%(grid_width **3))
    print ("Volume : %5.3f A^3"%(n_unoccupied*(grid_width**3)))
    return wrt,dat_wrt
def get_R_in_bsite(bsite_residues):
    res_R = []
    cutoff_list = []
    probe_radius = 1.4

    for atom in bsite_residues:
        atom_crd = atom.R
        element = get_element(atom.atm_name.strip())
        if not element in radius_dict:
            element = 'H'
        atom_radius = radius_dict[element]
        res_R.append(atom_crd.tolist())
        cutoff_list.append(atom_radius+probe_radius)

    res_R = torch.FloatTensor(res_R)
    cutoff_list = torch.FloatTensor(cutoff_list)
    return res_R, cutoff_list


def get_residues_in_bsite(protein, cntr_crd):
    bsite_residues = []
    for atom in protein.atom_s:
        dist = np.linalg.norm(atom.R-cntr_crd)
        if dist < 75:
            bsite_residues.append(atom.res_no)
    bsite_residues=list(set(bsite_residues))
    bsite_atom_s=[]
    for atom in protein.atom_s:
        if atom.res_no in bsite_residues:
            bsite_atom_s.append(atom)
    return bsite_atom_s


def get_protein_feature(protein,mat):
    cntr_crd=np.array([0.0,0.0,20.0]) # set z-coordinate of toggle residue C_alpha is 0
    bsite_residues = get_residues_in_bsite(protein, cntr_crd)
    volume_wrt,dat_wrt= check_volume_and_sasa(bsite_residues, cntr_crd,mat)
    return volume_wrt,dat_wrt





def main():
    if len (sys.argv)<2:
        print ("Enter 'python run.py -h' for print help message")
        print ("If you want to run demo,")
        print ("python run.py -p demo_input/demo_5zbq_trim.pdb -toggle 276 -tip_s 39,103,109,177,205,289,295")
        sys.exit()
    opt = argparse.ArgumentParser()

    opt.add_argument('-p','--pdb_fn',dest='pdb_fn',required=True,\
            help=' A query pdb file. ex) -p demo_5zbq.pdb')
    opt.add_argument('-toggle','--toggle_resno',dest='toggle_resno',required=True,type=int,\
            help=' Residue number of toggle switch residue. ex) -toggle 276')
    opt.add_argument('-tip_s','--tip_resno_s',dest='tip_resno_s',required=True,\
            help=' Tip residues of TM1~7.ex) -tip_s 39,103,109,177,205,289,295')
    opt.add_argument('-trim_Nterm','--trim_Nterm',dest='trim_Nterm',required=False,\
            type=int,default=35,\
            help=' Exclude residues, which have lower residue number than this variable,\
            for calculating center of protein. default :35')
    opt.add_argument('-trim_Cterm','--trim_Cterm',dest='trim_Cterm',required=False,\
            type=int,default=9999,\
            help=' Exclude residues, which have larger residue number than this variable,\
            for calculating center of protein. default :9999')
    opt.add_argument('-exclude_TM1_side_truncation','--exlcude_tm1',dest='exclude_tm1',required=False,\
            default=True,type=bool,\
            help=' Exclude TM1 for side truncation. When you think a TM1 is quite far from center of binding cavity\
            and or volume region between TM1,2,7 is separated from main binding cavity, you should turn this option\
            True ex. -exclude_TM1 False. default:True')
    #
    fn = opt.parse_args(); tag=fn.pdb_fn.split('/')[-1][:-4]
    tip_resno_s= [ int(x.strip()) for x in fn.tip_resno_s.split(',')]
    #
    prot=PDB(fn.pdb_fn)
    prot=trs_to_cntr(prot,fn.toggle_resno,fn.trim_Nterm,fn.trim_Cterm)
    cnt_s,helix_wrt=calc_helix_center(prot,tip_resno_s)
    tip_s=get_tips(cnt_s)
    mat=gen_mat()
    mat=truncate_side(mat,cnt_s,fn.exclude_tm1)
    mat,truncated_wrt=truncate_roof(mat,tip_s)
    #
    volume_wrt,dat_wrt=get_protein_feature(prot,mat)
    prot.write_pdb('%s_aligned_rec.pdb'%tag)
    with open('%s_helix.pdb'%tag,'wt')as fp:
        fp.writelines(helix_wrt)
    with open('%s_truncated.pdb'%tag,'wt')as fp:
        fp.writelines(truncated_wrt)
    with open('%s_volume.pdb'%tag,'wt')as fp:
        fp.writelines(volume_wrt)
    with open('%s_data.txt'%tag,'wt')as fp:
        fp.writelines(dat_wrt)
if __name__ =='__main__':
    main()
