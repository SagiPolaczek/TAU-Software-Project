from setuptools import setup, Extension

"""
A minimalist setup is shown.
"""


setup(name='capi_demo3',
      version='1.0',
      description='capi_demo 3 for sp class',
      ext_modules=[Extension('capi_demo3', sources=['capi3.c'])])
