from setuptools import find_packages, setup

setup(
    name="umichee",
    version="0.0.1",
    keywords=("pip", "umichee"),
    description="UMIche",
    long_description="Simulation tool UMIche",
    license="GNU GENERAL V3.0",

    url="https://github.com/cribbslab/umichee",
    author="Jianfeng Sun",
    # author_email="jianfeng.sun@ndorms.ox.ac.uk; Adam.Cribbs@ndorms.ox.ac.uk",

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
            'umichee=umichee.main:main',
        ],
    },
)