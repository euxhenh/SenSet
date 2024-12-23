import os

from setuptools import Command, find_packages, setup


class CleanCommand(Command):
    user_options = []  # type: ignore

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        os.system('rm -vrf ./build ./dist ./*.pyc ./*.tgz ./*.egg-info')


cmdclass = {'clean': CleanCommand}

setup(
    name='senset',
    packages=find_packages('.'),
    provides=['senset'],
    license='MIT',
    package_dir={'': '.'},
    version="0.0.1",
    cmdclass=cmdclass,
    author="Euxhen Hasanaj",
    author_email="ehasanaj@cs.cmu.edu",
    description=("Senescence Marker list in the Lung"),
    long_description=("Long Senescence Marker list in the Lung"),
    python_requires=">=3.10",
)
