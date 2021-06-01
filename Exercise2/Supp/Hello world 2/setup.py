from distutils.core import setup, Extension

"""
The old way of doing things, using distutils.
In addition, a minimalist setup is shown.
"""


setup(name='capi_demo2',
      version='1.0',
      description='capi_demo2 for sp class',
      ext_modules=[Extension('capi_demo2', sources=['capi2.c'])])
