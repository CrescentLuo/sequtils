import setuptools
 
with open('README.md', 'r') as fh:
    long_description = fh.read()
 
setuptools.setup(
    name='sequtils',
    version='0.0.1',
    author='Zheng Luo',
    author_email='crescentluozheng@gmail.com',
    description='Utils to high-throughput sequencing analysis',
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=setuptools.find_packages(),
    url='https://github.com/CrescentLuo/sequtils',
    classifiers=[
        'Programming Language :: Python :: 3',
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    keywords = 'sequencing',
    python_requires='>=3.6',
)