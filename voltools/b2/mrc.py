import mrcfile as mrc
import numpy as np


def save_mrc(data: np.ndarray, path: str, dtype=np.float32):
    f = mrc.new(path)

    # convert if needed
    if data.dtype != dtype:
        data = data.astype(dtype)

    # save
    f.set_data(data)

    # close
    f.close()


def read_mrc(path: str):
    f = mrc.open(path, permissive=True)
    data = f.data
    f.close()

    return data
