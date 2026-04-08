# CellBender Environment

Pattern: `hybrid`

Config variable: `CELLBENDER_ENV`

## Conda setup

```bash
source config/paths.sh
conda env create \
  -f Preprocessing/envs/cellbender/environment.yml \
  -p "${CELLBENDER_ENV}"
```

## Requirements setup

Create the CUDA-enabled base environment first, then install the Python layer:

```bash
source config/paths.sh
conda create -y -p "${CELLBENDER_ENV}" \
  -c conda-forge -c pytorch -c defaults \
  python=3.10 pip pytorch>=2.0 pytorch-cuda>=11.8
conda run -p "${CELLBENDER_ENV}" python -m pip install \
  -r Preprocessing/envs/cellbender/requirements.txt
```

## Notes

- Requires a GPU with enough VRAM for CellBender.
- Match the CUDA runtime to your cluster driver.
