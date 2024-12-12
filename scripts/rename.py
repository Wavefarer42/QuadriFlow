import pathlib
import shutil


def rename_files_in_folder(path_in: pathlib.Path, path_out: pathlib.Path):
    path_out.mkdir(parents=True, exist_ok=True)
    try:
        files_in = [it for it in sorted(path_in.iterdir()) if it.suffix in ['.ubs', '.png']]
        num_digits = len(str(len(files_in)))

        i = 0
        for path_file in files_in:
            filename_new = f"Z-{path_file.stem[:2]}{path_file.suffix}"
            path_new = path_out.joinpath(filename_new)
            shutil.copy(path_file, path_new)

            i += 1

        print(f"Renamed {num_digits} files successfully.")

    except Exception as e:
        print(f"An error occurred: {e}")


# Specify the folder path here
path_in = pathlib.Path('/Users/hannes/Projects/QuadriFlow/tests/resources/benchmark')
path_out = pathlib.Path('/Users/hannes/Projects/QuadriFlow/tests/resources/benchmark-out')
rename_files_in_folder(path_in, path_out)
