import os
import pathlib

cmd = "cmake-build-release/meshbound"

pathlib.Path("examples", exist_ok=True)
os.system(f"{cmd} --input tests/resources/sphere.ply --output examples/sphere-10k-ne-nb-na.obj --faces 10000")
os.system(f"{cmd} --input tests/resources/sphere.ply --output examples/sphere-10k-e-b-a.obj --faces 10000 --edges --adaptive")
os.system(f"{cmd} --input tests/resources/sphere.ply --output examples/sphere-5k-ne-nb-na.obj --faces 5000")
os.system(f"{cmd} --input tests/resources/sphere.ply --output examples/sphere-1k-ne-nb-na.obj --faces 1000")
os.system(f"{cmd} --input tests/resources/box.ply --output examples/box-10k-ne-nb-na.obj --faces 10000")
os.system(f"{cmd} --input tests/resources/box.ply --output examples/box-10k-e-b-a.obj --faces 10000 --edges --adaptive")
os.system(f"{cmd} --input tests/resources/box.ply --output examples/box-5k-ne-nb-na.obj --faces 5000")
os.system(f"{cmd} --input tests/resources/box.ply --output examples/box-1k-ne-nb-na.obj --faces 1000")
os.system(f"{cmd} --input tests/resources/box.ply --output examples/box-1k-e-nb-na.obj --faces 1000 --edges")
os.system(f"{cmd} --input tests/resources/fandisk.obj --output examples/fandisk-10k-ne-nb-na.obj --faces 10000")
os.system(f"{cmd} --input tests/resources/fandisk.obj --output examples/fandisk-10k-e-b-a.obj --faces 10000 --edges --adaptive")
os.system(f"{cmd} --input tests/resources/fandisk.obj --output examples/fandisk-5k-ne-nb-na.obj --faces 5000")
os.system(f"{cmd} --input tests/resources/fandisk.obj --output examples/fandisk-1k-ne-nb-na.obj --faces 1000")
os.system(f"{cmd} --input /Users/hannes/Data/datasets/botijo.ply --output examples/botijo-10k-ne-nb-na.obj --faces 10000")
os.system(f"{cmd} --input /Users/hannes/Data/datasets/botijo.ply --output examples/botijo-10k-e-b-a.obj --faces 10000 --edges --adaptive")
os.system(f"{cmd} --input /Users/hannes/Data/datasets/botijo.ply --output examples/botijo-5k-ne-nb-na.obj --faces 5000")
os.system(f"{cmd} --input /Users/hannes/Data/datasets/botijo.ply --output examples/botijo-1k-ne-nb-na.obj --faces 1000")
os.system(f"{cmd} --input /Users/hannes/Data/datasets/carter.ply --output examples/carter-10k-ne-nb-na.obj --faces 10000")
os.system(f"{cmd} --input /Users/hannes/Data/datasets/carter.ply --output examples/carter-10k-e-b-a.obj --faces 10000 --edges --adaptive")
os.system(f"{cmd} --input /Users/hannes/Data/datasets/carter.ply --output examples/carter-5k-ne-nb-na.obj --faces 5000")
os.system(f"{cmd} --input /Users/hannes/Data/datasets/carter.ply --output examples/carter-1k-ne-nb-na.obj --faces 1000")
