import argparse
import csv
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from tqdm import tqdm

def read_p_csv(path):
    p_lis = []
    with path.open('r') as f:
        reader = csv.reader(f)
        for row in reader:
            p_lis.append([round(float(p), 3) for p in row])
    return np.array(p_lis)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    paths = sorted(Path(args.input).glob('*.csv'))

    x_np = np.linspace(-10, 30, 401)
    y_np = np.linspace(-10, 10, 201)
    x, y = np.meshgrid(x_np, y_np)
    
    out_dir = Path(args.output).resolve()
    if not out_dir.exists():
        out_dir.mkdir()

    for i, path in enumerate(tqdm(paths)):
        arr = read_p_csv(path)
        plt.imshow(arr, extent=[np.min(x), np.max(x),
                                np.min(y), np.max(y)], vmin=-0.1, vmax=0.1)
        plt.savefig(str(out_dir.joinpath(f"{str(i).zfill(6)}.png")))
        plt.cla()

if __name__ == "__main__":
    main()
