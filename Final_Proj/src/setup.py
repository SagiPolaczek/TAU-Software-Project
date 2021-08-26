from distutils.core import setup, Extension

"""
The old way of doing things, using distutils.
In addition, a minimalist setup is shown.
"""


setup(name='myspkmeans',
      version='4.2',
      description='Our SAVAGE Program!',
      ext_modules=[Extension('myspkmeans', sources=['spkmeansmodule.c', 'spkmeans.c'])])
