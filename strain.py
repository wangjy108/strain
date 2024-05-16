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
#from util.Cluster import cluster


class main():
    def __init__(self, **args):
        try:
            self.db_name = args["input_sdf"]
        except Exception as e:
            self.db_name = None
        
        self.main_dir = os.getcwd()
        
        ## initial
        _dic_rotableBond = {0: 80, 1: 150} ## used to be 100, 300
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
        

        """
        try:
            self.rmsd_cutoff_initial_opt = float(args["rmsd_cutoff_initial_opt"])
        except Exception as e:
            self.rmsd_cutoff_initial_opt = 0.5
        
        """
        
        #print(f"gen conformer: {self.N_gen_conformer}")
        
        #self.HA_constrain = args["use_constrain"]
        #self.if_csv = args["if_csv"]  

        ##parameter for sample & cluster
        """
        try:
            self.k = args["cluster_n"]
        except Exception as e:
            self.k = None
        else:
            if not isinstance(self.k, int):
                self.k = None
        """

        
        try:
            self.distance_cutoff = args["rmsd_cutoff_cluster"]
        except Exception as e:
            self.distance_cutoff = 0.85
        else:
            if not isinstance(self.distance_cutoff, float):
                self.distance_cutoff = 0.85
        
        try:
            self.energy_cutoff = args["energy_cutoff"]
        except Exception as e:
            self.energy_cutoff = 0.5
        else:
            if not isinstance(self.energy_cutoff, float):
                self.energy_cutoff = 0.5
        
        try:
            self.sp_max_n = int(args["n_conformer_for_sp"])
        except Exception as e:
            self.sp_max_n = 1
               
        try:
            self.verbose = args["verbose"]
        except Exception as e:
            self.verbose = False

        try:
            self.if_double_constraint = args["double_constraint"]
        except Exception as e:
            self.if_double_constraint = False
            

        self.prefix = ".".join(self.db_name.split(".")[:-1])
      
    
    def step1_initial_opt(self):
        ## try without HeavyAtom constrain first
        #opt = []

        if self.charge != None:
            _initial_opt = sysopt(input_sdf=self.db_name, 
                                define_charge=self.charge,
                                HA_constrain=True).run()
        else:
            _initial_opt = sysopt(input_sdf=self.db_name,
                                  HA_constrain=True).run()
            self.charge = _initial_opt[0].GetProp("charge")
        
        """
        get_aligned_search = align(SearchMolObj=_initial_opt[0], RefMolObj=self.get_mol, method="crippen3D").run()
        if get_aligned_search:
            get_rmsd = RMSD(rdmolobj_mol=_initial_opt[0], rdmolobj_ref=self.get_mol, method="selfWhole").run()
        else:
            get_rmsd = None
            
        if get_rmsd > self.rmsd_cutoff_initial_opt:
            logging.info("Shift to HA constrained optimization")
            opt = sysopt(input_sdf=self.db_name,
                            HA_constrain=True, 
                            define_charge=self.charge).run()

        else:
            logging.info("Initial opt met rmsd requirement, no constrain performed")
            opt = _initial_opt
        """
        
        #if not opt:
        if not _initial_opt:
            logging.info("Failed at initial opt")
            return None
        
        cc = Chem.SDWriter(f"{self.prefix}.initial_opt.sdf")
        ## since this one the HA constrained conformation, should be ok if extra-align doesn't perform
        cc.write(_initial_opt[0])
        cc.close()

        return f"{self.prefix}.initial_opt.sdf"
    
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
    
    """
    def step3_cluster(self, input_sdf_name):

        #filtered_mol = cluster(input_sdf=input_sdf_name,
        #        cluster_n=self.k,
        #        rmsd_cutoff_cluster=self.distance_cutoff,
        #        if_write_cluster_mol=False).run()
        
        #filtered_file_name = [cc for cc in os.listdir() if ("FILTER_" in cc) and (cc.endswith(".sdf"))]
        if not filtered_mol:
            logging.info("Failed at cluster")
            return None
        
        cc = Chem.SDWriter("FILTER.sdf")
        for each in filtered_mol:
            cc.write(each)
        cc.close()

        return "FILTER.sdf"
    """
    
    def step3_opt(self, input_sdf_name):

        tracker = []
        
        opted = sysopt(input_sdf=input_sdf_name, 
                       rmsd_cutoff=self.distance_cutoff,
                      # save_n=self.sp_max_n,
                       HA_constrain=self.if_double_constraint,
                       define_charge=self.charge,
                       if_write_sdf=False).run()
        
        if not opted:
            logging.info("Failed at optimization")
            return None

        sorted_opt = sorted(opted, key=lambda x: float(x.GetProp("Energy_xtb")))
        i = 1
        saved = sorted_opt[:1]

        while i < len(sorted_opt):

            energy_diff = [abs(float(saved[ii].GetProp("Energy_xtb"))- float(sorted_opt[i].GetProp("Energy_xtb"))) * 627.51
                            for ii in range(len(saved))]

            if min(energy_diff) > self.energy_cutoff:
                saved.append(sorted_opt[i])

            i += 1
        
        return saved
    
    def step4_sp(self, mol_set):
        ## decide whether carry out sp calc

        if self.special_type:
                applied_basis = "6-311G**"
        else:
            if int(self.charge) < 0:
                applied_basis = "6-311G**"
            else:
                applied_basis = "6-31G**"

        
        mol_initial_opt = [cc for cc in Chem.SDMolSupplier(f"{self.prefix}.initial_opt.sdf", removeHs=False) if cc]
        xtb_initial_opt = float(mol_initial_opt[0].GetProp("Energy_xtb"))

        pose_with_energy_label_initial_opt, _ = syssp(input_rdmol_obj=mol_initial_opt, 
                                                charge_method="define", 
                                                define_charge=self.charge, 
                                                basis=applied_basis,
                                                if_solvent=True).run_pyscf()

        mol_initial_opt[0].SetProp("Energy_dft", f"{pose_with_energy_label_initial_opt[0][0]}")
        starter = Chem.SDWriter(f"{self.prefix}.initial_opt.sdf")
        starter.write(mol_initial_opt[0])
        starter.close()

        self.sp_max_n = min(len(mol_set), self.sp_max_n)

        flag = 0

        sp_dic = {}
        
        for ii, each_mol in enumerate(mol_set[:self.sp_max_n]):
            pose_with_energy_label, _ = syssp(input_rdmol_obj=[each_mol], 
                                                charge_method="define", 
                                                define_charge=self.charge, 
                                                basis=applied_basis,
                                                if_solvent=True).run_pyscf()

            each_mol.SetProp("Energy_dft", f"{pose_with_energy_label[0][0]}")
            ee_strain = (pose_with_energy_label_initial_opt[0][0] - pose_with_energy_label[0][0]) * 627.51

            sp_dic.setdefault(ii, [each_mol, ee_strain])
        
        get_GM = sorted(sp_dic.items(), key=lambda x: x[1][-1])[-1]
        
        #if get_GM[-1][-1] > 0:
        #    GM = mol_initial_opt[0]
        #    GM_energy = "~ 1"
        #else:
        #    GM = get_GM[-1][0]
        #    GM_energy = abs(get_GM[-1][-1])

        GM = get_GM[-1][0]
        GM_energy = get_GM[-1][-1]
        
        if self.sp_max_n > 1:
            more_cc = Chem.SDWriter(f"AppendixConf_top{self.sp_max_n}.sdf")
            for each in sorted(sp_dic.items(), key=lambda x: x[1][-1]):
                aligned_each = align(SearchMolObj=each[-1][0], RefMolObj=mol_initial_opt[0], method="crippen3D").run()
                more_cc.write(aligned_each)
            more_cc.close()
            logging.info(f"All conformations with DFT level energy labels were saved in [AppendixConf_top{self.sp_max_n}.sdf]")

        if self.verbose:
            more_and_more_cc = Chem.SDWriter(f"AppendixConf_rest.sdf")
            for mol in mol_set[self.sp_max_n:]:
                aligned_mol = align(SearchMolObj=mol, RefMolObj=mol_initial_opt[0], method="crippen3D").run() 
                more_and_more_cc.write(aligned_mol)
            more_and_more_cc.close()
            logging.info(f"More conformations without DFT energy labels were saved in [AppendixConf_rest.sdf]")
        
        
        if GM_energy < 0:
            write_GM_energy = "< 1"
            GM = mol_initial_opt[0]
        else:
            write_GM_energy = GM_energy
            ## perform align
            GM = align(SearchMolObj=GM, RefMolObj=mol_initial_opt[0], method="crippen3D").run() 

        cc = Chem.SDWriter("Global.sdf")
        cc.write(GM)
        cc.close()
        logging.info("Global minima conformation has been saved in [Global.sdf]")
            

        df = pd.DataFrame({"LigandName": [self.prefix],
                            "Applied Basis": [applied_basis],
                            "Strain Energy/kcal.mol-1": [f"{write_GM_energy}"]})


        df.to_csv("RESULT.STRAIN.csv", index=None)
        logging.info("Recorded strain energy saved in [RESULT.STRAIN.csv]")

        #os.system(f"rm -f OPT_*.sdf SAVE.sdf FILTER.sdf")

        return "RESULT.STRAIN.csv"
    
    def run(self):
        if not self.db_name:
            logging.info("check and run again")
            return

        work_dir = os.path.join(self.main_dir, f"{self.prefix}")
        if not os.path.exists(work_dir):
            os.mkdir(work_dir)
        os.chdir(work_dir)
        target = os.path.join(self.main_dir, self.db_name)
        os.system(f"mv {target} {work_dir}")

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
        
        #logging.info("Step3: Clustering ...")
        #self.filter_name = self.step3_cluster(self.sampled_file)
        #if not self.filter_name:
        #    logging.info("Terminated at Step3.")
        #    return
        
        logging.info("Step3: Optimization ...")
        self.opted_mol = self.step3_opt(self.sampled_file)
        if not self.opted_mol:
            logging.info("Terminated at Step4.")
        
        logging.info("Step4: Energy calculation ...")
    
        result_file_tag = self.step4_sp(self.opted_mol)
        if result_file_tag:
            logging.info("Finish strain calculation")
        
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
                        if not defined, 80 if rotable bond <=5, else 150')
    parser.add_argument('--define_charge', type=str, default=None,
                        help='define charge for input, use you sense that rdkit may fail to get corret bond order, default None')
    #parser.add_argument('--rmsd_cutoff_initial_opt', type=float,
    #                    help='define rmsd cutoff for initial opt, default is 0.5')
    #parser.add_argument('--not_save_csv', default=True, action='store_false', \
    #                    help="adding this option will not save final result in csv file")
    #parser.add_argument('--cluster_n', type=int, 
    #                    help='set static n for clustering after sampling')
    parser.add_argument('--rmsd_cutoff', type=float, default=0.85,
                        help='rmsd cutoff for reduce redundancy, default 0.85')
    parser.add_argument('--n_conformer_for_sp', type=int, default=1,
                        help='define the maximum conformers involved in final dft calculation, default 1')
    parser.add_argument('--energy_cutoff', type=float, default=0.5,
                        help='energy cutoff to save pose for single point calculation, unit in kcal/mol, default 0.5')
    parser.add_argument('--double_constraint', default=False, action='store_true',
                        help='if perform HeavyAtom constraint for post sampling optmization, default False')
    parser.add_argument('--verbose', default=False, action='store_true',
                        help='if save more conformation for visual check after single point energy calculation, default False')
    
    #parser.add_argument('--no_constrain', default=True, action='store_false', \
    #                    help="adding this option will turnoff constrain on heavy atoms when performing geometry optimization")
    args = parser.parse_args()

    main(input_sdf=args.input_sdf, 
         N_gen_conformer=args.N_gen_conformer, 
         define_charge=args.define_charge, 
         rmsd_cutoff_cluster=args.rmsd_cutoff,
         n_conformer_for_sp=args.n_conformer_for_sp,
         energy_cutoff=args.energy_cutoff,
         double_constraint=args.double_constraint,
         verbose=args.verbose).run()

    
##

