import pathlib
import shutil


def rename_files_in_folder(path_in: pathlib.Path, path_out: pathlib.Path):
    path_out.mkdir(parents=True, exist_ok=True)
    try:
        dir_cases: list[pathlib.Path] = [it for it in sorted(path_in.iterdir()) if it.is_dir()]

        for dir_case in dir_cases:
            name = dir_case.name
            path_file = [it for it in sorted(dir_case.iterdir())][-1]
            shutil.copy(path_file, path_out.joinpath(f"{name}{path_file.suffix}"))

    except Exception as e:
        print(f"An error occurred: {e}")


# Specify the folder path here
path_in = pathlib.Path('/Users/hannes/Projects/QuadriFlow/tests/out/benchmark')
path_out = pathlib.Path('/Users/hannes/Projects/QuadriFlow/tests/out/benchmark-package')
rename_files_in_folder(path_in, path_out)
