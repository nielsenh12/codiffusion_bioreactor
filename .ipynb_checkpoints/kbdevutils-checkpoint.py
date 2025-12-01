from __future__ import absolute_import

import sys
import json
import os
from os.path import exists
from pathlib import Path
from configparser import ConfigParser
from installed_clients.WorkspaceClient import Workspace
from kbbasemodules.basemodule import BaseModule
import logging

logger = logging.getLogger(__name__)
logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)

class KBDevUtils(BaseModule):
    def __init__(self,study_name,token_file=str(Path.home()) + '/.kbase/token',config_file=str(Path.home()) + '/.kbase/config',ws_version="prod",sdkhome=None,output_root=None):
        # Setting up filenames and checking for existence
        self.token_file = token_file
        if not exists(self.token_file):
            logger.critical("Token file not found! Use the util.set_token(<token>) to create the file!")
        self.config_file = config_file
        if not exists(self.config_file):
            logger.critical("You must create a config file in ~/.kbase/config before running this notebook. See instructions: https://docs.google.com/document/d/1fQ6iS_uaaZKbjWtw1MgzqilklttIibNO9XIIJWgxWKo/edit")
            sys.exit(1)
        # Loading config file
        confighash = self.read_config()
        # Setting up workspace URL
        self.ws_version = ws_version
        self.wsurl = self.get_wsurl_from_version(self.ws_version)
        # Setting up token
        token = self.read_token()
        # Setting up paths
        self.codebase = confighash.get("codebase","")
        # Running BaseModule init
        self.study_name = study_name
        # Getting callback url
        self.root_path = confighash.get("callback_path","")
        # Setting up SDK home
        self.sdkhome = sdkhome
        if not self.sdkhome:
            if self.ws_version == "prod":
                self.sdkhome = confighash.get("prod_sdk_home","prod")
            elif self.ws_version == "appdev":
                self.sdkhome = confighash.get("appdev_sdk_home","appdev")
        # Setting up output directory
        if output_root:
            self.output_root = output_root
        else:
            self.output_root = confighash.get("output_root","KBaseAnalysis/")
        if self.output_root[0] != "/":
            self.output_root = str(Path.home())+"/"+self.output_root
        self.output_dir = None
        self.data_dir = None
        self.set_study_root()
        print("Output files printed to:"+self.output_dir+" when using KBDevUtils.output_dir")
        # Setting callback directories
        self.callback_file = self.root_path+"/"+self.sdkhome+"/run_local/workdir/CallBack.txt"
        self.working_directory = self.root_path+"/"+self.sdkhome+"/run_local/workdir/tmp/"
        callback = self.read_callback()
        # Initializing BaseModule
        BaseModule.__init__(self,"KBDevUtils."+self.study_name,confighash,module_dir=self.codebase+"/chenry_utility_module/",working_dir=self.working_directory,token=token,clients={"Workspace":Workspace(self.wsurl, token=token)},callback=callback)
        self.version = "0.1.1.kbdu"
        self._msrecon = None
        self._kb_object_factory = None
    
    @property
    def msrecon(self):
        if self._msrecon is None:
            from ModelSEEDReconstruction.modelseedrecon import ModelSEEDRecon
            self._msrecon = ModelSEEDRecon(self.config,module_dir=self.codebase+"/KB-ModelSEEDReconstruction/",working_dir=self.working_dir,token=self.token,clients=self.clients,callback=self.callback_url)
        return self._msrecon
    
    @property
    def kb_object_factory(self):
        if self._kb_object_factory is None:
            from cobrakbase.core.kbase_object_factory import KBaseObjectFactory
            self._kb_object_factory = KBaseObjectFactory()
        return self._kb_object_factory

    def get_wsurl_from_version(self,version):
        if version == "prod":
            return "https://kbase.us/services/ws"
        elif version == "appdev":
            return "https://appdev.kbase.us/services/ws"
        elif version == "ci":
            return "https://ci.kbase.us/services/ws"
        else:
            logger.critical("Unknown workspace version: "+version)
            return "https://kbase.us/services/ws"

    def read_config(self):
        confighash = {}
        config = ConfigParser()
        config.read(self.config_file)
        for nameval in config.items("DevEnv"):
            confighash[nameval[0]] = nameval[1]
        return confighash

    def read_token(self):
        token = None
        with open(self.token_file, 'r') as fh:
            token = fh.read().strip()
        return token

    def read_callback(self):
        callback = None
        if exists(self.callback_file):
            with open(self.callback_file, 'r') as fh:
                callback = fh.read()
        return callback

    def check_kbase_dir(self):
        if not exists(str(Path.home()) + '/.kbase'):
            os.makedirs(str(Path.home()) + '/.kbase')
    
    def devutil_client(self):
        if "KBDevUtils" not in self.clients:
            from installed_clients.chenry_utility_moduleClient import chenry_utility_module
            self.clients["KBDevUtils"] = chenry_utility_module(self.callback_url,token=self.token)
        return self.clients["KBDevUtils"]
        
    def clear_sdk_dir(self):
        return self.devutil_client().run_command({"command":"clear"})
        
    def sdk_dir_perms(self):
        return self.devutil_client().run_command({"command":"perms"})
    
    def set_study_root(self,create=True):
        self.study_root = self.output_root+"/"+self.study_name+"/"
        if not exists(self.study_root):
            os.makedirs(self.study_root, exist_ok=True)
        self.data_dir = self.study_root + "datacache"
        os.makedirs(self.data_dir, exist_ok=True)
        self.output_dir = self.study_root + "nboutput"
        os.makedirs(self.output_dir, exist_ok=True)
    
    def set_token(self,token):
        self.check_kbase_dir()
        with open(self.token_file, 'w') as fh:
            fh.write(token)

    # Code for saving and loading data
    def save(self,name,data):
        filename = self.data_dir + "/" + name + ".json"
        dir = os.path.dirname(filename)
        os.makedirs(dir, exist_ok=True)
        with open(filename, 'w') as f:
            json.dump(data,f,indent=4,skipkeys=True)

    def load(self,name,default=None,kb_type=None):
        filename = self.data_dir + "/" + name + ".json"
        if not exists(filename):
            if default == None:
                logger.error("Requested data "+name+" doesn't exist at "+filename)
                raise(ValueError("Requested data "+name+" doesn't exist at "+filename))
            return default
        with open(filename, 'r') as f:
            data = json.load(f)
        if kb_type != None:
            data = self.kb_object_factory._build_object(kb_type, data, None, None)
        return data
        
    def list(self):
        files = os.listdir(self.data_dir)
        return [x.split(".")[0] for x in files if x.endswith(".json")]

    def exists(self,name):
        return exists(self.data_dir + "/" + name + ".json")

    # Code for printing object attributes and functions
    def object_attributes(self,obj,properties=True,functions=True):
        attributes = dir(obj)
        if properties:
            print("Properties:")
            properties = [attr for attr in attributes if not callable(getattr(obj, attr))]
            for property in properties:
                print(f"{property}")
        if functions:
            print("Functions:")
            functions = [attr for attr in attributes if callable(getattr(obj, attr))]
            for func in functions:
                print(f"{func}")