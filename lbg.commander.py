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
        parser.add_argument("--platform", help="local or lbg", default=False)
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
            self.main_dir = general["sdf_path"]
        except Exception as e:
            self.main_dir = None
        
        try:
            self.input_sdf_db = general["sdf_file"]
        except Exception as e:
            self.input_sdf_db = None
        
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
            self.rmsd_cutoff_cluster = strain["rmsd_cutoff"]
        except Exception as e:
            self.rmsd_cufoff_cluster = None
        
        if self.rmsd_cutoff_cluster:
            self._dic.setdefault("--rmsd_cutoff", self.rmsd_cutoff_cluster)
        
        #try:
        #    self.energy_cutoff = strain["energy_cutoff"]
        #except Exception as e:
        #    self.energy_cutoff = None
        
        #if self.energy_cutoff:
        #    self._dic.setdefault("--energy_cutoff", self.energy_cutoff)
        
        try:
            self.double_constraint = strain["double_constraint"]
        except Exception as e:
            self.double_constraint = False
        
        if self.double_constraint:
            self._dic.setdefault("--double_constraint", " ")
        
        try:
            self.verbose = strain["verbose"]
        except Exception as e:
            self.verbose = False
        
        if self.verbose:
            self._dic.setdefault("--verbose", " ")

        try:
            self.n_conformer_for_sp = strain["n_conformer_for_sp"]
        except Exception as e:
            self.n_conformer_for_sp = None
        
        if self.n_conformer_for_sp:
            self._dic.setdefault("--n_conformer_for_sp", self.n_conformer_for_sp)
        

        self.root_dir = os.getcwd()
        self.run_type = args.run_type
        self.platform = args.platform


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
    
    def read_sdf(self, db_name, save_path):
        rdmol_obj_dic = {}

        with open(db_name, "r+") as f:
            sdf_content = [ff for ff in f.readlines()]
        
        end_idx = [-1] + [idx for idx, line in enumerate(sdf_content) if "$$$$" in line][:-1]

        for ii, idx in enumerate(end_idx):
            try:
                get_mol = sdf_content[end_idx[ii]+1:end_idx[ii+1]+1]
            except Exception as e:
                get_mol = sdf_content[end_idx[ii]+1:]
            
            get_real_name = str(get_mol[0].strip())

            try:
                rdmol_obj_dic[get_real_name]
            except Exception as e:
                rdmol_obj_dic.setdefault(get_real_name, [])
            
            rdmol_obj_dic[get_real_name].append(get_mol)
        
        if not os.path.exists(save_path):
            os.mkdir(save_path)
        
        for kk, vv in rdmol_obj_dic.items():
            save_file = f"{os.path.join(save_path, kk)}.sdf"
            with open(save_file, "w+") as cc:
                for idx, each in enumerate(vv):
                    for line in each:
                        cc.write(line)
                 
        return [kk for kk in rdmol_obj_dic.keys()]

    def prepare(self):
        logging.info("Run prepare ....")

        if not (self.main_dir or self.input_sdf_db):
            logging.info("No available path/sdf to run, please define sdf_path or sdf_file")
            return
        
        if self.main_dir and (not os.path.exists(self.main_dir)):
            logging.info("Wrong path to run, please re-define sdf_path")
            return
        
        if self.input_sdf_db:
            self.main_dir = os.path.join(os.getcwd(), "STRAIN")
            _all = self.read_sdf(self.input_sdf_db, self.main_dir)
            os.chdir(self.main_dir)
        
        if self.main_dir:
            os.chdir(self.main_dir)
            _all = [sdf for sdf in os.listdir("./") if ".sdf" in sdf]

        if not _all:
            logging.info("Nothing to do, abort")
            return 
            
        cmd_py = os.path.join(self.root_dir, 'strain.py')
        
        if not os.path.isfile(cmd_py):
            logging.info(f"No strain.py, please put it at {self.root_dir}")
            return 
        
        if self.platform == "lbg":
            n_group = math.ceil(len(_all) / self.N_in_group)
        else:
            n_group = 1
            self.N_in_group = len(_all)

        
        _work_dir = [os.path.join(self.main_dir, str(i)) for i in range(n_group)]
        

        for ii, each_dir in enumerate(_work_dir):
            if not os.path.exists(each_dir):
                os.mkdir(each_dir)

            os.chdir(each_dir)
            
            cmd_cp = f"cp {cmd_py} ./"
            (_, _) = subprocess.getstatusoutput(cmd_cp)

            logging.info(f"Preparing at {each_dir}")

            get_set = _all[ii*self.N_in_group:(ii+1)*self.N_in_group]

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
                target_file = os.path.join(each_dir, each)
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
                                
                cmd_line = f"python strain.py --input_sdf {each} "

                if cmd_dic:
                    for kk, vv in cmd_dic.items():
                        cmd_line += f"{kk} {vv} "
                
                cmd_line += "\n"
                
                bash_cmd.append(cmd_line)

                if os.path.isfile(target_file):
                    _flag += 1

            if _flag == len(get_set):
                self.setup_local_special(work_dir=each_dir, cmd_list=bash_cmd)
            
            os.chdir(self.root_dir)
        
        logging.info("Finish prepare")

        return _work_dir

    def submit_lbg(self):
        tracker = {}
        if not self.lbg_project:
            logging.info("Please define lbg_project in config")
            return 
        
        if not self.lbg_image:
            logging.info("Please define lbg_image in config")
            return
        
        work_dir = self.prepare()

        logging.info("Submit job(s) to bohrium ....")
        
        for idx, _dir in enumerate(work_dir):
            os.chdir(_dir)
            self.setup_lbg(work_dir=_dir)

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
                    tracker.setdefault(_dir, JOB_ID)
                else:
                    logging.info(f"Submit failed at {_dir}, check in folder for more information")
            else:
                logging.info(f"Submit failed at {_dir}, check in folder for more information")

            os.chdir(self.root_dir)
        
        if tracker:
            df_track = pd.DataFrame({"local_path": [kk for kk in tracker.keys()],
                                 "JOB_ID": [vv for vv in tracker.values()]})
        
            df_track.to_csv(os.path.join(self.main_dir,"submission_track_list"), index=None)
        else:
            logging.info("Submit failed")
            return 

        logging.info("Submit finish")
        return 
    
    def submit_local(self):

        work_dir = self.prepare()

        logging.info("Run local strain calculation ....")

        tracker = {}
        for idx, _dir in enumerate(work_dir):
            os.chdir(_dir)

            logging.info(f"Working at {_dir}")

            submission_cmd = "bash run.sh"
            (status, out) = subprocess.getstatusoutput(submission_cmd)  

            tag = [line.strip().endswith("Done") for line in out.split("\n")]

            if tag:
                tracker.setdefault(_dir, "Done")
            else:
                tracker.setdefault(_dir, "Error")

            os.chdir(self.root_dir)
        
        if tracker:
            df_track = pd.DataFrame({"local_path": [kk for kk in tracker.keys()],
                                 "Status": [vv for vv in tracker.values()]})
        
            df_track.to_csv(os.path.join(self.main_dir,"submission_track_list"), index=None)
        else:
            logging.info("local run failed")
            return 

        logging.info("Local run finish")
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
        logging.info("Strain results for all calculated mols saved in [FinalStrain.csv] ")
        return 

    def col_result(self):
        if not self.main_dir:
            _dir = os.path.join(self.root_dir, "STRAIN")
            self.main_dir = _dir
        else:
            _dir = self.main_dir

        os.chdir(_dir)

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
            
            if "Status" in df.columns:
                if row["Status"] == "Done":
                    transfer += [os.path.join(row["local_path"], dd) for dd in os.listdir(row["local_path"]) if os.path.isdir(dd)]
                    should_delete.append(row["local_path"])
                else:
                    should_not_delet.append(row["local_path"])
            else:
                cmd = f"lbg job download {row['JOB_ID']}"
                (_, _) = subprocess.getstatusoutput(cmd)
                if os.path.isdir(str(row['JOB_ID'])):
                    (_, _) = subprocess.getstatusoutput(f"mv {row['JOB_ID']}/* ./")
                    (_, _) = subprocess.getstatusoutput(f"rm -rf {row['JOB_ID']}")
                    transfer += [os.path.join(row["local_path"], dd) for dd in os.listdir(row["local_path"]) if os.path.isdir(dd)]
                    should_delete.append(row["local_path"])
                else:
                    should_not_delet.append(row["local_path"])
            
            os.chdir(_dir)
        
        if transfer:
            target = os.path.join(_dir, "RUN_DONE")
            if not os.path.exists(target):
                os.mkdir(target)
            for each in transfer:
                cmd_mv = f"mv {each} {target}"
                os.system(cmd_mv)
            
            for pp in should_delete:
                shutil.rmtree(pp)
            
        if should_not_delet:
            error = os.path.join(_dir, "ERROR")
            if not os.path.exists(error):
                os.mkdir(os.path.join(_dir, "ERROR"))
            
            for _pp in should_not_delet:
                cmd_mv2 = f"mv {_pp}/*.sdf {error}"
                os.system(cmd_mv2)

        if os.path.exists("RUN_DONE"):
            logging.info(f"Finish collect")
            self.summarize(_path=os.path.join(_dir, "RUN_DONE"))
            os.system("rm -f submission_track_list")
        else:
            logging.info("Nothing collected, abort")
            return
        
        return
    
    def run(self):
        _dic_platform = {"lbg": self.submit_lbg,
                         "local": self.submit_local}

        if "submit" in self.run_type:
            try:
                _dic_platform[self.platform]
            except Exception as e:
                logging.info("Wrong platform option, choose from [local, lbg]")
                return

            _dic_platform[self.platform]()

        elif "col_result" in self.run_type:
            self.col_result()
            #os.system("rm -f submission_track_list")
        else:
            logging.info(f"Not available run type: {self.run_type}")
            logging.info("Choose from [submit, col_result]")
        
        return


if __name__ == "__main__":
    main().run()
    
    


    