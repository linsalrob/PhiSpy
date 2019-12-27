from distutils.core import setup, Extension
import io
import os
import modules.helper_functions

# this is taken from https://www.jeffknupp.com/blog/2013/08/16/open-sourcing-a-python-project-the-right-way/ and
# allows README.txt and CHANGES.txt to automatically be used by PyPI. Create README.txt using http://pandoc.org/:
# pandoc -t plain -o README.txt README.md


def read(*filenames, **kwargs):
    encoding = kwargs.get('encoding', 'utf-8')
    sep = kwargs.get('sep', '\n')
    buf = []
    for filename in filenames:
        if os.path.exists(filename):
            with io.open(filename, encoding=encoding) as f:
                buf.append(f.read())
    return sep.join(buf)

long_description = read('README.txt', 'CHANGES.txt')




def main():
    setup(
        name="PhiSpy",
        version=modules.helper_functions.get_version(),
        description="Prophage finder using multiple metrics",
        long_description=long_description,
        author="Rob Edwards",
        platforms='any',
        keywords="phage prophage bioinformatics microbiology bacteria genome genomics",
        author_email="raedwards@gmail.com",
        url='https://github.com/linsalrob/PhiSpy',
        license='The MIT License (MIT)',
        scripts=['PhiSpy.py'],
        zip_safe=True,
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
        entry_points = {
            'PhiSpy': [
                'eggsecutable = PhiSpy.py:main',
            ]
        }

        )

if __name__ == "__main__":
    main()
