import subprocess
from pathlib import Path

mzmine_path = Path(r"C:\Program Files\mzmine\mzmine_console.exe")


def run_mzmine_batch(batchfile):
    batch_path = Path(batchfile)

    if not mzmine_path.is_file():
        print(f"Error: mzmine_console.exe not found at {mzmine_path}")
        return 1

    if not batch_path.is_file():
        print(f"Error: Batch file not found at {batch_path}")
        return 1

    command = [str(mzmine_path), "-b", str(batch_path)]
    print(f"Running command: {' '.join(f'\"{arg}\"' if ' ' in arg else arg for arg in command)}")

    try:
        result = subprocess.run(command, check=False)
        print(f"Finished file {batchfile} with return code {result.returncode}")
        return result.returncode
    except OSError as exc:
        print(f"Error launching MZmine: {exc}")
        return 1




if __name__ == "__main__":
    batchfile = r"\\marty-nas.cm.utexas.edu\MartyLab\Group Share\Annika\DDA_Compare\mzmine\Stellar_Pos_batch_no_na.mzbatch"
    run_mzmine_batch(batchfile)
