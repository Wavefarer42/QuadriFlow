import itertools
import pathlib
import shutil
import subprocess

cmd = pathlib.Path("cmake-build-release/meshbound")

if __name__ == '__main__':
    dir_input = pathlib.Path("tests/resources/benchmark")
    dir_output = pathlib.Path("tests/out/benchmark")

    shutil.rmtree(dir_output, ignore_errors=True)
    dir_output.mkdir(parents=True, exist_ok=True)

    assert dir_input.exists(), "Input directory does not exist"
    assert cmd.exists(), "Meshbound executable does not exist"

    faces = [100, 1000, 10000]
    preserve_edges = [False, True]

    cases = list(itertools.product(faces, preserve_edges))

    for path_in in dir_input.iterdir():
        if path_in.suffix == ".ubs":
            for faces, preserve_edges in cases:
                arguments = [
                    cmd,
                    f"--input", path_in,
                    "--output", dir_output.joinpath(f"{path_in.stem}-{faces}-{'e' if preserve_edges else 'ne'}.ply"),
                    f"--faces", faces,
                    f"--edges" if preserve_edges else None
                ]

                subprocess.run([str(it) for it in arguments if it])
