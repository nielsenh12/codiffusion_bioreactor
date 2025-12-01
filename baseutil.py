import platform
import sys
import sys
import json
from json import dump
import os
import re
from os.path import exists
from pathlib import Path
import logging
import shutil

print("python version " + platform.python_version())

from configparser import ConfigParser

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

config = ConfigParser()
if not exists(str(Path.home()) + '/.kbase/config'):    
    if exists("/scratch/shared/code/sharedconfig.cfg"):
        shutil.copyfile("/scratch/shared/code/sharedconfig.cfg",str(Path.home()) + '/.kbase/config')
    else:
        print("You much create a config file in ~/.kbase/config before running this notebook. See instructions: https://docs.google.com/document/d/1fQ6iS_uaaZKbjWtw1MgzqilklttIibNO9XIIJWgxWKo/edit")
        sys.exit(1)
config.read(str(Path.home()) + '/.kbase/config')
paths = config.get("DevEnv","syspaths").split(";")
codebase = config.get("DevEnv","codebase",fallback="")
for i,filepath in enumerate(paths):
    if filepath[0:1] != "/":
        paths[i] = codebase+"/"+filepath
sys.path = paths + sys.path
display(sys.path)
from kbdevutils import KBDevUtils

class BaseUtil(KBDevUtils):
    def __init__(self,name):
        KBDevUtils.__init__(self,name,output_root=os.path.dirname(os.path.realpath(__file__)))
        
