import os
global xp
global device
global map_coordinates

if 'PYTOM_GPU' in os.environ.keys() and str(os.environ['PYTOM_GPU']) != '-1':
    try:
        import cupy as xp
        import cupy.typing as xptyping
        ID = os.environ['PYTOM_GPU'].split(',')
        xp.cuda.Device(int(ID[0])).use()
        from cupyx.scipy.ndimage import map_coordinates
        device = f'gpu:{ID[0]}'

    except Exception as e:
        print(e)
        import numpy as xp
        import numpy.typing as xptyping
        from scipy.ndimage import map_coordinates

        device = 'cpu'
else:
    os.system("export PYTOM_GPU=0")
    import numpy as xp
    import numpy.typing as xptyping
    from scipy.ndimage import map_coordinates

    device = 'cpu'

def initialize_gpu(id):
    global xp, device
    import cupy as xp
    xp.cuda.Device(int(id)).use()
    from cupyx.scipy.ndimage import map_coordinates
    device = f'gpu:{id}'
