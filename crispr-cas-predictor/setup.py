from setuptools import setup, find_packages

setup(
    name='crispr-cas-predictor',
    version='0.1.0',
    author='Your Name',
    author_email='your.email@example.com',
    description='A tool for predicting CRISPR-Cas proteins from protein sequences.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/yourusername/crispr-cas-predictor',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    install_requires=[
        'biopython',
        'numpy',
        'scikit-learn',
        'hmmlearn',
        'pandas',
        'matplotlib',
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)