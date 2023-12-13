import os
import json
import math
import subprocess
import argparse
import configparser
import logging
import shutil
import pandas as pd

logging.basicConfig(format='[%(levelname)s] %(message)s', level=logging.INFO)

class main():
    def __init__(self):
        parser = argparse.ArgumentParser(description="STRAIN")

        parser.add_argument("--config", help="path of config file", default=False)
        parser.add_argument("--run_type", help="submit or col_result", default=False)
        args = parser.parse_args()

        config = configparser.ConfigParser()
        config.read(args.config)

        self._dic = {}

        lbg_config = config["lbg_config"]

        try:
            self.lbg_project = lbg_config["lbg_project"]
        except Exception as e:
            self.lbg_project = None
        
        try:
            self.lbg_machine_type = lbg_config["lbg_machine_type"]
        except Exception as e:
            self.lbg_machine_type = "c16_m128_cpu"
        if not self.lbg_machine_type:
            self.lbg_machine_type = "c16_m128_cpu"
        
        try:
            self.lbg_image=lbg_config["lbg_image"]
        except Exception as e:
            self.lbg_image = None
        
        try:
            self.N_in_group = lbg_config.getint("N_in_group")
        except Exception as e:
            self.N_in_group = 50
        else:
            if not self.N_in_group:
                self.N_in_group = 50
        
        general = config["general"]

        try:
            self.main_dir = general["work_dir"]
        except Exception as e:
            self.main_dir = None
        
        try:
            self.extra_charge = general["extra_charge"]
        except Exception as e:
            self.extra_charge = None
        else:
            if not os.path.isfile(self.extra_charge):
                self.extra_charge = None
        
        strain = config["strain"]

        try:
            self.N_gen_conformer = strain["N_gen_conformer"]
        except Exception as e:
            self.N_gen_conformer = None
        
        if self.N_gen_conformer:
            self._dic.setdefault("--N_gen_conformer", self.N_gen_conformer)
        
        try:
            self.rmsd_cufoff_initial_opt = strain["rmsd_cutoff_initial_opt"]
        except Exception as e:
            self.rmsd_cufoff_initial_opt = None

        if self.rmsd_cufoff_initial_opt:
            self._dic.setdefault("--rmsd_cutoff_initial_opt", self.rmsd_cufoff_initial_opt)
        
        try:
            self.cluster_n = strain["cluster_n"]
        except Exception as e:
            self.cluster_n = None
        
        if self.cluster_n:
            self._dic.setdefault("--cluster_n", self.cluster_n)
        
        try:
            self.rmsd_cutoff_cluster = strain["rmsd_cutoff_cluster"]
        except Exception as e:
            self.rmsd_cufoff_cluster = None
        
        if self.rmsd_cutoff_cluster:
            self._dic.setdefault("--rmsd_cutoff", self.rmsd_cutoff_cluster)
        
        try:
            self.n_conformer_for_sp = strain["n_conformer_for_sp"]
        except Exception as e:
            self.n_conformer_for_sp = None
        
        if self.n_conformer_for_sp:
            self._dic.setdefault("--n_conformer_for_sp", self.n_conformer_for_sp)
        

        self.root_dir = os.getcwd()
        self.run_type = args.run_type


    def setup_lbg(self, work_dir):
        lbg_json = {
            "job_name": "STRAIN",
            "command": "bash run.sh",
            "log_file": "tmp_log",
            "backward_files": [],
            "program_id": f"{self.lbg_project}",
            "platform": "ali",
            "job_group_id": "",
            "disk_size": 128,
            "machine_type": f"{self.lbg_machine_type}",
            "job_type": "container",
            "image_name": f"{self.lbg_image}"
        }
        
        file_name = os.path.join(work_dir, "input.json")
        with open(file_name, "w+") as cc:
            json.dump(lbg_json, cc, indent=4)
        return

    
    def setup_local_special(self, work_dir, cmd_list):
        file_name = os.path.join(work_dir, "run.sh")
    
        with open(file_name, "w+") as cc:
            cc.write("#!/bin/sh \n")
            cc.write("\n")
            for each_line in cmd_list:
                cc.write(each_line)
            cc.write("\n")
        return

    def submit(self):
        if not self.main_dir:
            logging.info("No available path to run, please define work_dir in config")
            return 
        
        if not os.path.exists(self.main_dir):
            logging.info("Wrong path to run, please re-define work_dir in config")
            return
        
        os.chdir(self.main_dir)

        _all = [sdf for sdf in os.listdir("./") if ".sdf" in sdf]

        if not _all:
            logging.info(f"Nothing to run in {self.main_dir}, terminate")
            return
        
        if not self.lbg_project:
            logging.info("Please define lbg_project in config")
            return 
        
        if not self.lbg_image:
            logging.info("Please define lbg_image in config")
            return 


        cmd_py = os.path.join(self.root_dir, 'sample_pip_strain.py')
        
        if not os.path.isfile(cmd_py):
            logging.info(f"No sample_pip_strain.py, please put it at {self.root_dir}")
            return 
        
        logging.info("Run submit ....")
        
        n_group = math.ceil(len(_all) / self.N_in_group)

        tracker = {}

        for i in range(n_group):
            _work_dir = os.path.join(self.main_dir, str(i))
            if not os.path.exists(_work_dir):
                os.mkdir(_work_dir)
            os.chdir(_work_dir)
            logging.info(f"Working with subset {i}")
            get_set = _all[i*self.N_in_group:(i+1)*self.N_in_group]

            _flag = 0

            bash_cmd = []

            _dic_charge = {}

            if self.extra_charge:
                self.extra_charge = os.path.join(self.root_dir, self.extra_charge)
                df_charge = pd.read_csv(self.extra_charge)
                for idx, row in df_charge.iterrows():
                    _dic_charge.setdefault(row["LigLabel"], row["charge"])

            for each in get_set:
                original_file = os.path.join(self.main_dir, each)
                target_file = os.path.join(_work_dir, each)
                cmd = f"mv {original_file} {target_file}"
                (_, _) = subprocess.getstatusoutput(cmd)

                cmd_dic = {}
                
                try:
                    _dic_charge[each.split(".")[0]]
                except Exception as e:
                    cmd_dic.update(self._dic)
                else:
                    cmd_dic.update(self._dic)
                    cmd_dic.setdefault("--define_charge", _dic_charge[each.split(".")[0]])
                                
                cmd_line = f"python sample_pip_strain.py --input_sdf {each} "

                if cmd_dic:
                    for kk, vv in cmd_dic.items():
                        cmd_line += f"{kk} {vv} "
                
                cmd_line += "\n"
                
                bash_cmd.append(cmd_line)

                if os.path.isfile(target_file):
                    _flag += 1
            
            if _flag == len(get_set):
                #setup_lbg(projectID=projectID, _path=_work_dir)
                self.setup_lbg(work_dir=_work_dir)
                self.setup_local_special(work_dir=_work_dir, cmd_list=bash_cmd)
                
                cmd_cp = f"cp {cmd_py} ./"
                (_, _) = subprocess.getstatusoutput(cmd_cp)

                submission_cmd = "lbg job submit -i input.json -p ./ > TRACK_JOB"
                (_, _) = subprocess.getstatusoutput(submission_cmd)
            
            if os.path.isfile("TRACK_JOB") and os.path.getsize("TRACK_JOB"):
                with open("TRACK_JOB", "r+") as f:
                    content = [ll.strip() for ll in f.readlines() if ll]
                
                if content:
                    try:
                        JOB_ID = content[-1].split()[-1]
                    except Exception as e:
                        JOB_ID = None
                else:
                    JOB_ID = None
                
                if JOB_ID:
                    tracker.setdefault(_work_dir, JOB_ID)
                else:
                    logging.info(f"Submit failed at {_work_dir}, check in folder for more information")
            else:
                logging.info(f"Submit failed at {_work_dir}, check in folder for more information")

            os.chdir(self.main_dir)
        
        if tracker:
            df_track = pd.DataFrame({"local_path": [kk for kk in tracker.keys()],
                                 "JOB_ID": [vv for vv in tracker.values()]})
        
            df_track.to_csv(os.path.join(self.main_dir,"submission_track_list"), index=None)
        else:
            logging.info("Submit failed")
            return 

        logging.info("Finish submit")
        return 
    
    def get_path(self, root_path_prefix):
        if not isinstance(root_path_prefix, list):
            root_path_prefix = [root_path_prefix]
        new_root_path = []
        for each in root_path_prefix:
            new_root_path += [os.path.join(each, cc) for cc in os.listdir(each) if os.path.isdir(os.path.join(each, cc))]
        
        #print(new_root_path)

        if not new_root_path:
            return root_path_prefix
        else:
            return self.get_path(new_root_path)
    
    def summarize(self, _path):
        all_path = self.get_path(_path)
        avail = []
        nonavail = []
        for pp in all_path:
            #work_dir = os.path.join(_path, pp)
            csv = [cc for cc in os.listdir(pp) if cc.endswith("csv")]
            if csv:
                df = pd.read_csv(os.path.join(pp, csv[0]))
                #df["LigandName"] = df["mol_label"].apply(lambda x: "_".join(x.split("_")[:-1]))
                #df["conformer_ID"] = df["mol_label"].apply(lambda x: x.split("_")[-1])
                #df = df.drop(columns=["mol_label"])[["LigandName", "conformer_ID", "Strain_ene(kcal/mol)"]]
                avail.append(df)
            else:
                get_ligand_name = pp.split("/")[-1]
                nonavail.append(get_ligand_name)
        
        if avail:
            _df_all = pd.concat(avail)
            if nonavail:
                _df_append = pd.DataFrame(columns=_df_all.columns, index=[ii for ii in range(len(nonavail))])
                _df_append["LigandName"] = nonavail

                _df = pd.concat([_df_all, _df_append])
            else:
                _df = _df_all
        else:
            logging.info("Nothing to collect, can not process strain calculation for all input molecules")
            return
        
        _df.to_csv(os.path.join(self.main_dir, "FinalStrain.csv"), index=None)
        logging.info(f"Strain results for all calculated mols saved in [FinalStrain.csv] at {self.main_dir}")
        return 

    def col_result(self):
        os.chdir(self.main_dir)

        if not os.path.isfile("submission_track_list"):
            logging.info("No tracking information available, abort")
            return 

        logging.info("Collecting results ....")

        df = pd.read_csv("submission_track_list", header=0)

        should_not_delet = []
        transfer = []
        should_delete = []

        for idx, row in df.iterrows():
            os.chdir(row["local_path"])
            cmd = f"lbg job download {row['JOB_ID']}"
            (_, _) = subprocess.getstatusoutput(cmd)
            if os.path.isdir(str(row['JOB_ID'])):
                (_, _) = subprocess.getstatusoutput(f"mv {row['JOB_ID']}/* ./")
                (_, _) = subprocess.getstatusoutput(f"rm -rf {row['JOB_ID']}")
                transfer += [os.path.join(row["local_path"], dd) for dd in os.listdir(row["local_path"]) if os.path.isdir(dd)]
                should_delete.append(row["local_path"])
            else:
                should_not_delet.append(row["local_path"])
            
            os.chdir(self.main_dir)
        
        if transfer:
            target = os.path.join(self.main_dir, "RUN_DONE")
            if not os.path.exists(target):
                os.mkdir(target)
            for each in transfer:
                cmd_mv = f"mv {each} {target}"
                os.system(cmd_mv)
            
            for pp in should_delete:
                shutil.rmtree(pp)
            
        if should_not_delet:
            error = os.path.join(self.main_dir, "ERROR")
            if not os.path.exists(error):
                os.mkdir(os.path.join(self.main_dir, "ERROR"))
            
            for _pp in should_not_delet:
                cmd_mv2 = f"mv {_pp}/*.sdf {error}"
                os.system(cmd_mv2)

        if os.path.exists("RUN_DONE"):
            logging.info("Finish collect")
            self.summarize(_path=os.path.join(self.main_dir, "RUN_DONE"))
            os.system("rm -f submission_track_list")
        else:
            logging.info("Nothing collected, abort")
            return
        
        return
    
    def run(self):
        if "submit" in self.run_type:
            self.submit()
        elif "col_result" in self.run_type:
            self.col_result()
            #os.system("rm -f submission_track_list")
        else:
            logging.info(f"Not available run type: {self.run_type}")
        
        return


if __name__ == "__main__":
    main().run()
    
    


    