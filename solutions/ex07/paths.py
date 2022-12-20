from typing import Optional
import pathlib

from global_utils import results_dir, ROOT

INPUT_DIR = ROOT / "solutions" / "ex07"


def zero_dir(section: int, phase: str, step: str, method: Optional[str] = None) -> pathlib.Path:
    if method is None:
        return INPUT_DIR / str(section) / phase / step
    else:
        return INPUT_DIR / str(section) / phase / step / method


def measures_dir(section: int, phase: str, step: str, method: Optional[str] = None) -> pathlib.Path:
    """Directory where simulation results for the given phase and step will be stored
    """
    if method is None:
        return results_dir("07") / str(section) / phase / step
    else:
        return results_dir("07") / str(section) / phase / step / method
