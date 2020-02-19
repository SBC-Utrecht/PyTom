import os

if 'PYTOM_GPU' in os.environ.keys() and os.environ['PYTOM_GPU']:
    try:
        import cupy as xp
        device='gpu'
        print('GPU code activated')

    except:
        import numpy as xp
        device='cpu'
else:
    os.system("export PYTOM_GPU=0")
    import numpy as xp
    device='cpu'

