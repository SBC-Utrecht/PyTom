from pytom.gpu.initialize import xp
def volumesSameSize(volume1: xp.typing.NDArray, volume2: xp.typing.NDArray) -> bool:
    return volume1.shape == volume2.shape
