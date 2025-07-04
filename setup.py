import setuptools
from distutils.core import Extension
import os

# Read the markdown files for the long description
def read(*filenames, **kwargs):
    encoding = kwargs.get('encoding', 'utf-8')
    sep = kwargs.get('sep', '\n')
    buf = []
    for filename in filenames:
        if os.path.exists(filename):
            with open(filename, "r", encoding=encoding) as f:
                buf.append(f.read())
    return sep.join(buf)

long_description = read('README.md', 'CHANGES.md')

def get_version():
    with open("VERSION", 'r') as f:
        v = f.readline().strip()
        return v

def get_requirements():
    reqs = []
    with open('requirements.txt', 'r') as f:
        for l in f:
            reqs.append(l.strip())
    return reqs

def main():
    setuptools.setup(
        name="PhiSpy",
        version=get_version(),
        description="Prophage finder using multiple metrics",
        long_description=long_description,
        long_description_content_type="text/markdown",
        author="Rob Edwards",
        platforms='any',
        keywords="phage prophage bioinformatics microbiology bacteria genome genomics",
        author_email="raedwards@gmail.com",
        url='https://github.com/linsalrob/PhiSpy',
        license='The MIT License (MIT)',
        scripts=['PhiSpy.py', 'scripts/make_training_sets.py', 'scripts/plot_stats.py', 'scripts/compare_predictions_to_phages.py', 'scripts/mark_prophage_features.py'],
        packages=setuptools.find_packages(),
        ext_modules=[Extension("PhiSpyRepeatFinder", sources=["src/repeatFinder.cpp"], language='c++')],
        include_package_data=True,
        classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: MIT License',
            'Natural Language :: English',
            'Operating System :: MacOS :: MacOS X',
            'Operating System :: POSIX :: Linux',
            'Operating System :: Unix',
            'Programming Language :: Python :: 3.0',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
        install_requires = [
            'biopython >= 1.74',
            'numpy >= 1.16.0',
            'scikit-learn >= 0.21.3',
            'bcbio-gff >= 0.6.6'
        ],
        entry_points={
            "console_scripts": ["PhiSpy.py = PhiSpyModules.main:run", "phispy = PhiSpyModules.main:run"]
        }

    )

if __name__ == "__main__":
    main()
