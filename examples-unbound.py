import os
import pathlib

cmd = "cmake-build-release/meshbound"

pathlib.Path("examples").mkdir(exist_ok=True)
os.system(f"{cmd} --input /Users/hannes/Data/Unbound/models/Demon.ubs --output examples/demon-10k-ne-nb-na.obj --faces 10000")
os.system(f"{cmd} --input /Users/hannes/Data/Unbound/models/Greenbot_001.ubs --output examples/greenbot-10k-ne-nb-na.obj --faces 10000")
os.system(f"{cmd} --input /Users/hannes/Data/Unbound/models/model_head.ubs --output examples/model_head-10k-ne-nb-na.obj --faces 10000")
os.system(f"{cmd} --input /Users/hannes/Data/Unbound/models/Pinball.ubs --output examples/pinball-10k-ne-nb-na.obj --faces 10000")
os.system(f"{cmd} --input /Users/hannes/Data/Unbound/models/valheim.ubs --output examples/valheim-10k-ne-nb-na.obj --faces 10000")
os.system(f"{cmd} --input /Users/hannes/Data/Unbound/models/wizards.ubs --output examples/wizards-10k-ne-nb-na.obj --faces 10000")
