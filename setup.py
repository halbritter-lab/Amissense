import os
from setuptools import setup, find_packages

# Load version from version.py
version = {}
with open(os.path.join('Amissense', 'version.py')) as f:
    exec(f.read(), version)

setup(
    name="Amissense",
    version=version['__version__'],
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "biopython>=1.84",
        "requests>=2.28.0",
        "pandas>=2.0.0",
        "numpy>=1.23.0",
        "seaborn>=0.12.0",
        "matplotlib>=3.6.0",
        "plotly>=5.11.0",
        "kaleido>=0.2.1",
    ],
    entry_points={
        "console_scripts": [
            "amissense=Amissense.cli:main",
        ],
    },
    author="Claudia Abad Baucells",
    author_email="claudia.abad-baucells@charite.de",
    description="Amissense: A tool for processing AlphaMissense data and generating visualizations",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/halbritter-lab/Amissense",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',
    extras_require={
        "dev": [
            "pytest",
            "black",
            "flake8",
        ],
    },
)
