#!/usr/bin/env python
# coding: utf-8

"""
// author: Wang (Max) Jiayue 
// email: wangjy108@outlook.com
"""

import argparse
import logging
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import os
import sys
import pandas as pd
import numpy as np
import shutil
from util.ConfGenbyMM import ConfGen
from util.OptbySQM import System as sysopt
from util.SPcalc import System as syssp
from util.Align import Align as align
from util.CalRMSD import RMSD
from util.ConfRelaxbySQM import System as MDsample
from util.Cluster import cluster


class main():
    def __init__(self, **args):
        try:
            self.db_name = args["input_sdf"]
        except Exception as e:
            self.db_name = None
        
        self.main_dir = os.getcwd()
        
        ## initial
        _dic_rotableBond = {0: 100, 1: 300}
        #_dic_HA_index = {}
        try:
            self.get_mol = [mm for mm in Chem.SDMolSupplier(self.db_name) if mm][0]
        except Exception as e:
            self.get_mol = None
            rotable_bond_index = 0
        else:
            rotable_bond = rdMolDescriptors.CalcNumRotatableBonds(self.get_mol)
            _def_func = lambda x: 0 if max(5, x) == 5 else 1
            rotable_bond_index = _def_func(rotable_bond)
        
        self.special_type = False

        if self.get_mol:
            atom_type = set([atom.GetSymbol() for atom in self.get_mol.GetAtoms()])
            if "I" in atom_type:
                self.special_type = True
        
        try:
            self.N_gen_conformer = args["N_gen_conformer"]
        except Exception as e:
            self.N_gen_conformer = _dic_rotableBond[rotable_bond_index]    
        else:
            try:
                self.N_gen_conformer = int(self.N_gen_conformer)
            except Exception as e:
                self.N_gen_conformer = _dic_rotableBond[rotable_bond_index]
        
        try:
            self.charge = args["define_charge"]
        except Exception as e:
            self.charge = None
        else:
            try:
                self.charge = int(self.charge)
            except Exception as e:
                self.charge = None
        
        try:
            self.rmsd_cutoff_initial_opt = float(args["rmsd_cutoff_initial_opt"])
        except Exception as e:
            self.rmsd_cutoff_initial_opt = 0.5
        
        #print(f"gen conformer: {self.N_gen_conformer}")
        
        #self.HA_constrain = args["use_constrain"]
        #self.if_csv = args["if_csv"]  

        ##parameter for sample & cluster
        try:
            self.k = args["cluster_n"]
        except Exception as e:
            self.k = None
        else:
            if not isinstance(self.k, int):
                self.k = None

        try:
            self.distance_cutoff = args["rmsd_cutoff_cluster"]
        except Exception as e:
            self.distance_cutoff = 0.85
        else:
            if not isinstance(self.distance_cutoff, float):
                self.distance_cutoff = 0.85
        
        try:
            self.sp_max_n = args["n_conformer_for_sp"]
        except Exception as e:
            self.sp_max_n = 1
      
    
    def step1_initial_opt(self):
        ## try without HeavyAtom constrain first
        opt = []
        if self.charge != None:
            _initial_opt = sysopt(input_sdf=self.db_name,
                                 charge_method="define", 
                                 define_charge=self.charge).run()
        else:
            _initial_opt = sysopt(input_sdf=self.db_name).run()
            self.charge = int(_initial_opt[0].GetProp("charge"))
        
        get_aligned_search = align(SearchMolObj=_initial_opt[0], RefMolObj=self.get_mol, method="crippen3D").run()
        if get_aligned_search:
            get_rmsd = RMSD(rdmolobj_mol=_initial_opt[0], rdmolobj_ref=self.get_mol, method="selfWhole").run()
        else:
            get_rmsd = None
            
        if get_rmsd > self.rmsd_cutoff_initial_opt:
            logging.info("Shift to HA constrained optimization")
            opt = sysopt(input_sdf=self.db_name,
                            HA_constrain=True, 
                            charge_method="define", 
                            define_charge=self.charge).run()
        else:
            logging.info("Initial opt met rmsd requirement, no constrain performed")
            opt = _initial_opt
        
        if not opt:
            logging.info("Failed at initial opt")
            return None
        
        cc = Chem.SDWriter(f"{self.db_name.split('.')[0]}.initial_opt.sdf")
        cc.write(opt[0])
        cc.close()

        return f"{self.db_name.split('.')[0]}.initial_opt.sdf"
    
    def step2_sampling(self, input_file_name):
        sampled_file_name = None
        MDsample(input_sdf=input_file_name, 
                 save_frame=self.N_gen_conformer, 
                 define_charge=self.charge).run()
        if os.path.isfile("SAVE.sdf") and os.path.getsize("SAVE.sdf"):
            logging.info("Samping success")
            sampled_file_name = "SAVE.sdf"
        
        else:
            logging.info("Failed at sampling")

        return sampled_file_name
    
    def step3_cluster(self, input_sdf_name):

        filtered_mol = cluster(input_sdf=input_sdf_name,
                cluster_n=self.k,
                rmsd_cutoff_cluster=self.distance_cutoff,
                if_write_cluster_mol=False).run()
        
        #filtered_file_name = [cc for cc in os.listdir() if ("FILTER_" in cc) and (cc.endswith(".sdf"))]
        if not filtered_mol:
            logging.info("Failed at cluster")
            return None
        
        cc = Chem.SDWriter("FILTER.sdf")
        for each in filtered_mol:
            cc.write(each)
        cc.close()

        return "FILTER.sdf"
    
    def step4_opt(self, input_sdf_name):

        tracker = []
        
        opted = sysopt(input_sdf=input_sdf_name, 
                       rmsd_cutoff=self.distance_cutoff,
                       save_n=self.sp_max_n,
                       HA_constrain=False, 
                       define_charge=self.charge,
                       if_write_sdf=False).run()
        
        if not opted:
            logging.info("Failed at optimization")
            return None

        for idx, mm in enumerate(opted):
            cc = Chem.SDWriter(f"OPT_{idx}.sdf")
            cc.write(mm)
            cc.close()
            tracker.append(f"OPT_{idx}.sdf")
        
        return tracker
    
    def step5_sp(self, tracker):
        ## decide whether carry out sp calc

        if self.special_type:
                applied_basis = "6-311G*"
        else:
            if self.charge < 0:
                applied_basis = "6-31+G*"
            else:
                applied_basis = "6-311G*"

        
        mol_initial_opt = [cc for cc in Chem.SDMolSupplier(f"{self.db_name.split('.')[0]}.initial_opt.sdf", removeHs=False) if cc]
        xtb_initial_opt = float(mol_initial_opt[0].GetProp("Energy_xtb"))

        opted_mol = {}

        flag = 0

        for idx, name in enumerate(tracker):
            get_mol =  [cc for cc in Chem.SDMolSupplier(name, removeHs=False) if cc]
            get_energy = (float(get_mol[0].GetProp("Energy_xtb")) - xtb_initial_opt) * 627.51
            if get_energy < flag:
                flag = get_energy
            opted_mol.setdefault(name, get_mol)
        
        if flag == 0:
            ## no need for sp calc
            GM = mol_initial_opt
            GM_energy = "~ 1"
            write_basis = None
        else:
            sp_dic = {}
            pose_with_energy_label_initial_opt, _ = syssp(input_rdmol_obj=mol_initial_opt, 
                                                 charge_method="define", 
                                                 define_charge=self.charge, 
                                                 basis=applied_basis).run_pyscf()
            
            mol_initial_opt[0].SetProp("Energy_dft", f"{pose_with_energy_label_initial_opt[0][0]}")

            for kk, vv in opted_mol.items():
                pose_with_energy_label, _ = syssp(input_rdmol_obj=vv, 
                                                 charge_method="define", 
                                                 define_charge=self.charge, 
                                                 basis=applied_basis).run_pyscf()
                
                vv[0].SetProp("Energy_dft", f"{pose_with_energy_label[0][0]}")
                ee_strain = (pose_with_energy_label[0][0] - pose_with_energy_label_initial_opt[0][0]) * 627.51

                sp_dic.setdefault(kk, [vv, ee_strain])
            
            get_GM = sorted(sp_dic.items(), key=lambda x: x[1][-1])[0]
            
            if get_GM[-1][-1] > 0:
                GM = mol_initial_opt
                GM_energy = "~ 1"
            else:
                GM = get_GM[-1][0]
                GM_energy = abs(get_GM[-1][-1])
            
        cc = Chem.SDWriter(f"Stable_{self.db_name.split('.')[0]}.sdf")
        cc.write(GM[0])
        cc.close()

        df = pd.DataFrame({"LigandName": [self.db_name.split('.')[0]],
                            "Applied Basis": [applied_basis],
                            "Strain Energy/kcal.mol-1": [f"{GM_energy}"]})

        df.to_csv("RESULT.STRAIN.csv", index=None)

        #os.system(f"rm -f OPT_*.sdf SAVE.sdf FILTER.sdf")

        return "RESULT.STRAIN.csv"
    
    def run(self):
        if not self.db_name:
            logging.info("check and run again")
            return
        
        self.prefix = ".".join(self.db_name.split(".")[:-1])
        work_dir = os.path.join(self.main_dir, f"{self.prefix}")
        if not os.path.exists(work_dir):
            os.mkdir(work_dir)
        os.chdir(work_dir)
        os.system(f"mv {self.main_dir}/{self.db_name} ./")

        ## step1
        logging.info("Step1: Performing initial opt ....")

        self.initial_opt = self.step1_initial_opt()

        if not self.initial_opt:
            logging.info("Terminated at Step1.")
            return 
        
        ## step2
        logging.info("Step2: Geometrical Sampling ....")
        self.sampled_file = self.step2_sampling(self.initial_opt)
        if not self.sampled_file:
            logging.info("Terminated at Step2.")
            return 
        
        logging.info("Step3: Clustering ...")
        self.filter_name = self.step3_cluster(self.sampled_file)
        if not self.filter_name:
            logging.info("Terminated at Step3.")
            return
        
        logging.info("Step4: Optimization ...")
        self.opted_name_list = self.step4_opt(self.filter_name)
        if not self.opted_name_list:
            logging.info("Terminated at Step4.")
        
        logging.info("Step5: Energy calculation ...")
        result_file_tag = self.step5_sp(self.opted_name_list)
        if result_file_tag:
            logging.info("Finish strain calculation")
            logging.info("Stable pose save in Stable_*.sdf")
            logging.info("Others in *.csv")
        
        logging.info("Done")
        
        os.system(f"rm -f OPT_*.sdf SAVE.sdf FILTER.sdf")

        os.chdir(self.main_dir)




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='sample naive strain energy workflow, \
                                    $WORKDIR has same name with input sdf, \
                                    save energy labeled sdf in [**_withEneTag.sdf], \
                                    and stable pose in [opt_***.sdf]')
    parser.add_argument('--input_sdf', type=str, required=True, 
                        help='input sdf file, docking pose(s) for single mol')
    parser.add_argument('--N_gen_conformer', type=int, 
                        help='available for [denovo and sample] mode, define N conformers in sampling stage, \
                        if not defined, 100 if rotable bond <=5, else 300')
    parser.add_argument('--define_charge', type=str, default=None,
                        help='define charge for input, use you sense that rdkit may fail to get corret bond order, default None')
    parser.add_argument('--rmsd_cutoff_initial_opt', type=float,
                        help='define rmsd cutoff for initial opt, default is 0.5')
    #parser.add_argument('--not_save_csv', default=True, action='store_false', \
    #                    help="adding this option will not save final result in csv file")
    parser.add_argument('--cluster_n', type=int, 
                        help='set static n for clustering after sampling')
    parser.add_argument('--rmsd_cutoff', type=float,
                        help='rmsd cutoff for reduce redundancy, default 0.85')
    parser.add_argument('--n_conformer_for_sp', type=int,
                        help='define the maximum conformers involved in final dft calculation, default 1')
    #parser.add_argument('--no_constrain', default=True, action='store_false', \
    #                    help="adding this option will turnoff constrain on heavy atoms when performing geometry optimization")
    args = parser.parse_args()

    main(input_sdf=args.input_sdf, \
         N_gen_conformer=args.N_gen_conformer, \
         define_charge=args.define_charge, \
         rmsd_cutoff_initial_opt=args.rmsd_cutoff_initial_opt,\
         clsuter_n=args.cluster_n, \
         rmsd_cutoff_cluster=args.rmsd_cutoff, \
         n_conformer_for_sp=args.n_conformer_for_sp).run()




    
##

