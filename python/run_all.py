# runs everything in order: python python/run_all.py
import runpy, pathlib

here = pathlib.Path(__file__).parent
for step in ("analysis.py", "figures.py", "make_kml.py"):
    print(f"--- {step}")
    runpy.run_path(str(here / step), run_name="__main__")
