# let guys could run Goodvibes by python -m goodvibes <args>
# Copied from __main__.py in pip
from __future__ import absolute_import

import os
import sys

# If we are running from a wheel, add the wheel to sys.path
# This allows the usage python pip-*.whl/pip install pip-*.whl

if __package__ == '':
    path = os.path.dirname(os.path.dirname(__file__))
    sys.path.insert(0, path)

from readgv import readgv  # noqa

if __name__ == '__main__':
    readgv.main()
#    sys.exit(vasputils.main()) 