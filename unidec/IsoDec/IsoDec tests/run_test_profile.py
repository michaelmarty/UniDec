import os

import torch

from unidec.IsoDec.engine import IsoDecEngine
from unidec.IsoDec.train_profile import get_model_input_settings
from unidec.IsoDec.encoding_profile import data_dirs


def load_test_data(dirs, phaseres, phaselock, mode):
    """Load test data for a trained profile model.

    Parameters
    ----------
    dirs : iterable of path-like
        Directories containing ``training_data_profile.npz`` and
        ``test_data_profile.npz``. These are the same directory inputs accepted
        by :func:`unidec.IsoDec.train_profile.run_training`.
    phaseres : int
        Phase resolution used by the trained model.
    phaselock : int
        Phase lock used by the trained model.
    mode : str
        Mode used by the trained model. Must be one of ``"centroid"``,
        ``"profile"``, or ``"both"``.

    Returns
    -------
    IsoDecEngine
        Engine with the test data loaded.
    """
    eng = IsoDecEngine(phaseres=int(phaseres))
    original_dir = os.getcwd()
    try:
        eng.create_merged_dataloader(
            [os.fspath(directory) for directory in dirs],
            "profile",
            noise_percent=0,
            batchsize=1024,
            double_percent=0.0,
            onedrop_percent=0.0,
            harmonic_percent=0,
            equilize=False,
            ematkey=get_model_input_settings(phaseres, phaselock, mode)[0],
            centroidkey=get_model_input_settings(phaseres, phaselock, mode)[1],
            dataset="test"
        )
    finally:
        os.chdir(original_dir)
    return eng


def get_model_accuracy(model_input, dirs, model_dir=None, model_override=None, topdir=None, save_test_data_flag=True):
    """Evaluate one trained profile model against the training script's test data.

    Parameters
    ----------
    model_input : tuple
        ``(phaseres, phaselock, mode)``, where mode is ``"centroid"``,
        ``"profile"``, or ``"both"``.
    dirs : iterable of path-like
        Directories containing ``training_data_profile.npz`` and
        ``test_data_profile.npz``. These are the same directory inputs accepted
        by :func:`unidec.IsoDec.train_profile.run_training`.
    model_dir : path-like, optional
        Directory containing the trained model. Defaults to IsoDec's
        ``modelparams`` directory.

    Returns
    -------
    float
        Fraction of correct predictions in the range 0 to 1.
    """
    try:
        phaseres, phaselock, mode = model_input
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "model_input must contain (phaseres, phaselock, mode)"
        ) from exc

    ematkey, centroidkey, model_name = get_model_input_settings(phaseres, phaselock, mode)
    if model_override is not None:
        model_name = "phase_model_" + model_override + ".pth"

    if model_dir is None:
        model_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "modelparams")
    model_path = os.path.abspath(os.path.join(os.fspath(model_dir), model_name))
    if not os.path.isfile(model_path):
        raise FileNotFoundError(f"Model not found: {model_path}")

    if isinstance(dirs, (str, os.PathLike)):
        dirs = [dirs]

    if save_test_data_flag:
        # Check if save file already exists
        outfile = os.path.join(topdir if topdir is not None else os.getcwd(),
                               f"test_dataloader_{phaseres}_{phaselock}_{mode}.pth")
        if os.path.isfile(outfile):
            eng = load_test_data_from_file(phaseres, phaselock, mode, topdir=topdir)
        else:
            eng = load_test_data(dirs, phaseres, phaselock, mode)
            save_test_data(eng, phaseres, phaselock, mode, overwrite=True, topdir=topdir)
    else:
        eng = load_test_data(dirs, phaseres, phaselock, mode)

    eng.phasemodel.setup_training(forcenew=False)
    eng.phasemodel.load_model(model_path)
    _, accuracy = eng.phasemodel.evaluate_model(eng.test_dataloader)
    return accuracy


def save_test_data(eng, phaseres, phaselock, mode, overwrite=False, topdir=None):
    """Save the test data to a file.

    Parameters
    ----------
    eng: IsoDecEngine
        Engine with the test data loaded.
    phaseres : int
        Phase resolution.
    phaselock : int
        Phase lock.
    mode : str
        Mode of operation.

    Returns
    -------
    None
    """
    # Check for if file exists
    fileout = f"test_dataloader_{phaseres}_{phaselock}_{mode}.pth"
    if os.path.isfile(fileout) and not overwrite:
        print(f"File already exists: {fileout}. Use overwrite=True to overwrite.")
        return
    elif os.path.isfile(fileout) and overwrite:
        print(f"Overwriting existing file: {fileout}")
    elif not os.path.isfile(fileout):
        print(f"Saving test data to new file: {fileout}")

    outfile = os.path.join(topdir if topdir is not None else os.getcwd(), fileout)
    torch.save(eng.test_dataloader, outfile)


def load_test_data_from_file(phaseres, phaselock, mode, topdir=None):
    """Load the test data from a file.

    Parameters
    ----------
    phaseres : int
        Phase resolution.
    phaselock : int
        Phase lock.
    mode : str
        Mode of operation.

    Returns
    -------
    IsoDecEngine
        Engine with the test data loaded.
    """
    eng = IsoDecEngine(phaseres=int(phaseres))
    outfile = os.path.join(topdir if topdir is not None else os.getcwd(),
                           f"test_dataloader_{phaseres}_{phaselock}_{mode}.pth")
    test_dataloader = torch.load(outfile, weights_only=False)
    eng.test_dataloader = test_dataloader
    return eng


def test_grid(phaseres, dirs, topdir=None):
    print("Running Test Grid for Phase Resolution:", phaseres)
    if phaseres > 1:
        model_inputs = [
            (phaseres, 0, "profile"),
            (phaseres, 1, "profile"),
            (phaseres, 0, "centroid"),
            (phaseres, 1, "centroid"),
            (phaseres, 0, "both"),
            (phaseres, 1, "both"), ]

        model_overrides = [
            "p_{}_0".format(phaseres),
            "p_{}_1".format(phaseres),
            "c_{}_0".format(phaseres),
            "c_{}_1".format(phaseres),
            "b_{}_0".format(phaseres),
            "b_{}_1".format(phaseres),
            "{}".format(phaseres)
        ]
    else:
        model_inputs = [
            (phaseres, 0, "profile"),
            (phaseres, 0, "centroid"),
            (phaseres, 0, "both"), ]

        model_overrides = [
            "p_{}_0".format(phaseres),
            "c_{}_0".format(phaseres),
            "b_{}_0".format(phaseres),
            "{}".format(phaseres)
        ]

    results = []
    for inp in model_inputs:
        for override in model_overrides:
            try:
                accuracy = get_model_accuracy(inp, dirs, model_override=override, topdir=topdir)
                print(f"Model input: {inp}, Override: {override}, Accuracy: {accuracy:.4f}")
                results.append((inp, override, accuracy))
            except FileNotFoundError as e:
                print(f"Model input: {inp}, Override: {override}, Error: {e}")

    # Write to output text file in output directory
    output_file = os.path.join(topdir if topdir is not None else os.getcwd(), f"model_accuracy_results_{phaseres}.tsv")
    with open(output_file, "w") as f:
        f.write("Model Input\tOverride\tAccuracy\n")
        for inp, override, accuracy in results:
            f.write(f"{inp}\t{override}\t{accuracy:.4f}\n")

if __name__ == "__main__":
    topdirectory = "/groups/mtmarty/data"
    topdirectory = "Z:\\Group Share\\JGP\\"
    testdir = "Z:\\Group Share\\JGP\\TestData\\"

    # dirs = [os.path.join(topdirectory, d) for d in data_dirs]
    #
    # test_grid(4, dirs=dirs, topdir=topdirectory)
    #
    # exit()
    # Example usage
    model_input = (8, 1, "profile")

    # dirs = ["path/to/training_data", "path/to/test_data"]
    dirs = ["PXD069439"]
    accuracy = get_model_accuracy(model_input, dirs, model_override="b_8_1", topdir=topdirectory)
    print(f"Model accuracy: {accuracy:.4f}")
