#!/bin/python
# -*- coding:utf-8 -*-
import os
import numpy as np
import shutil

class VASP(object):
    def __init__(self):
        self.lattice=np.zeros((3,3))
        self.Selective_infomation=False
        self.Direct_mode=True

    def xyz_read(self):
        #print('Now reading vasp structures.')
        if os.path.exists("POSCAR"):
            poscar=open("POSCAR",'r')
        elif os.path.exists("CONTCAR"):
            poscar=open("CONTCAR",'r')
        else:
            raise IOError('CONTCAR OR POSCAR does not exist!')

        self.title=poscar.readline().rstrip('\r\n').rstrip('\n')
        self.scaling_factor=float(poscar.readline())
        for i in range(3):
            self.lattice[i]=np.array([float(j) for j in poscar.readline().split()])
        #self.lattice*=self.scaling_factor
        self.element_list=[j for j in poscar.readline().split()]
        try:
            self.element_amount=[int(j) for j in poscar.readline().split()]
        except ValueError:
            raise ValueError('VASP 5.x POSCAR is needed!')
        line_tmp=poscar.readline()
        if line_tmp.strip().upper().startswith("S"):
            self.Selective_infomation=True 
            line_tmp=poscar.readline()  
        else:# no atoms fixed
            self.Selective_infomation=False 
        if line_tmp.strip().upper().startswith("D"):
            self.Direct_mode=True
        elif line_tmp.strip().upper().startswith("C"):
            self.Direct_mode=False
        else:
            raise ValueError("POSCAR format is not correct!")
        total_atom=sum(self.element_amount)
        self.atomic_position=np.zeros((total_atom,3))
        self.Selective_TF=[]
        if self.Selective_infomation == True:   
            for i in range(total_atom):
                line_tmp=poscar.readline()
                self.atomic_position[i]=np.array([float(j) for j in line_tmp.split()[0:3]])
                self.Selective_TF.append([j for j in line_tmp.split()[3:]]) 
        else:
            for i in range(total_atom):
                line_tmp=poscar.readline()
                self.atomic_position[i]=np.array([float(j) for j in line_tmp.split()[0:3]])         
        if self.Direct_mode == True:
            self.atomic_position=np.dot(self.atomic_position,self.lattice)
            self.atomic_position*=self.scaling_factor
        poscar.close()
        self.lattice*=self.scaling_factor
        return (np.matmul(self.atomic_position,np.linalg.inv(self.lattice)),self.lattice,self.element_list,self.element_amount)

if __name__ == "__main__":
    poscar=VASP()
    atoms,lattice,element_list,element_amount=poscar.xyz_read()
    with open("sphere.ini",'w') as writer:
        writer.write("[geometry]\nlattice_vectors: \n")
        writer.write(" %f %f %f\n" %(lattice[0][0],lattice[0][1],lattice[0][2]))
        writer.write(" %f %f %f\n" %(lattice[1][0],lattice[1][1],lattice[1][2]))
        writer.write(" %f %f %f\n\nbasis:\n" %(lattice[2][0],lattice[2][1],lattice[2][2]))
        count = 0 
        for index,_ in enumerate(element_list):
            for i in range(element_amount[index]):
                writer.write(" %s %f %f %f\n" %(element_list[index],atoms[count][0],atoms[count][1],atoms[count][2]))
                count+=1
        writer.write("basis_coordsys: lattice\n\nshift_vector: 0.0 0.0 0.0\nshift_vector_coordsys: lattice\n\n")
        writer.write("[sphere: 1]\n\nradius: 10\n")
    
