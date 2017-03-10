import os
from setuptools import setup, find_packages

#
# For packing data see https://danielsokolowski.blogspot.ru/2012/08/setuptools-includepackagedata-option.html
# $ python3 setup.py sdist
#


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(
    name="pystella",
    version="0.1",
    author="Petr Baklanov",
    author_email="baklanovp@gmail.com",
    description="To work with supernovae models computed by STELLA.",
    license="MIT",
    keywords="supernova, light curves",
    url="https://github.com/baklanovp/pystella",
    packages=find_packages(),
    # packages=['pystella', 'pystella.fit', 'pystella.model', 'pystella.rf', 'pystella.util'],
    # data_files=[('', ['data/*'])],
    # package_data={'': ['data/']},
    include_package_data=True,
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Programming Language :: Python :: 3.5",
        "License :: OSI Approved :: MIT License",
    ],
)
