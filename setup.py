from setuptools import find_packages, setup

setup(
    name="umiche",
    version="0.0.1",
    keywords=["conda", "umiche"],
    description="UMIche",
    long_description="Simulation tool UMIche",
    license="GNU GENERAL V3.0",

    url="https://github.com/cribbslab/umiche, https://github.com/2003100127",
    author="Jianfeng Sun",
    author_email="jianfeng.sun@ndorms.ox.ac.uk, adam.cribbs@ndorms.ox.ac.uk",

    packages=find_packages(),
    include_package_data=True,
    # package_data={},
    platforms="any",
    python_requires=">3.8",
    install_requires=[
        "click",
        "pandas",
        "numpy",
        "scipy",
        "biopython",
        "pyyaml",
        "pyfiglet==0.8.post1",
    ],
    entry_points={
        'console_scripts': [
            'umiche=umiche.main:main',
        ],
    },
)