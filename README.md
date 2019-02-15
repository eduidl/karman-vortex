# Karman Vortex

## Requirement

- cmake
- Python 3.6
- pipenv

## Run

```sh
mkdir build && cd build
cmake .. && make && ./karman
cd ../
pipenv sync
pipenv run python3 csv2png.py --input data --output images
ffmpeg -r 120 -i ../images/%06d.png -vcodec libx264 -pix_fmt yuv420p -r 120 result.mp4
```
