import os
import sys


currDir  = os.path.dirname(os.path.abspath(__file__))
projRoot = os.path.abspath(os.path.join(currDir, ".."))
srcPath  = os.path.join(projRoot, "src")

if (srcPath not in sys.path):

    sys.path.insert(0, srcPath)