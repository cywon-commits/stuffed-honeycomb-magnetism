import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from analyze_crossing import analyze
from refine_crossings import summarize
for j in [.35,.425]:
 d=summarize(j);assert d['global_ground_identified'];analyze(d['winner']['name'])
