"""Run fixture tests, then execute every analysis-notebook cell in the local kernel."""
import ast
import os
import sys
import unittest
from pathlib import Path

import nbformat
from nbclient import NotebookClient

ROOT = Path(__file__).resolve().parents[1]
os.chdir(ROOT)
sys.path.insert(0, str(ROOT))
suite = unittest.defaultTestLoader.discover(str(ROOT / 'tests'))
if not unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful():
    raise SystemExit(1)
path = ROOT / 'dog_ocd_analysis.ipynb'
notebook = nbformat.read(path, as_version=4)
for cell in notebook.cells:
    if cell.cell_type == 'code':
        ast.parse(cell.source)
pilot = '--pilot' in sys.argv
if pilot:
    notebook.cells = notebook.cells[:28]
def report(cell, cell_index, **kwargs):
    print(f'Executing cell {cell_index}: {cell.source.splitlines()[0][:100]}', flush=True)
client = NotebookClient(notebook, timeout=1800, kernel_name='python3',
                        resources={'metadata': {'path': str(ROOT)}}, on_cell_execute=report)
try:
    client.execute()
finally:
    destination = ROOT / 'notebook/results/ipynb/toga/pilot.ipynb' if pilot else path
    destination.parent.mkdir(parents=True, exist_ok=True)
    nbformat.write(notebook, destination)
print('Executed', destination, flush=True)
