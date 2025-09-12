from setuptools import Extension, setup

module = Extension("symnmf_mod", sources=['symnmfmodule.c', 'symnmf.c'])
setup(name='symnmf_mod',
     version='1.0',
     description='Python wrapper for custom C extension',
     ext_modules=[module])