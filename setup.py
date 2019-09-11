from setuptools import setup

setup(
        name="hic_hiker",
        version="1.0.0",
        description="3D-DNA scaffolds refinement",
        long_description=open('README.md').read(),
        long_description_content_type="text/markdown",
        url="https://github.com/ryought/hic_hiker",
        author="ryought",
        author_email="ryonakabayashi@gmail.com",
        license="MIT",
        packages=[
            "hic_hiker",
            ],
        install_requires=[
            "numpy",
            "matplotlib",
            "scikit-learn",
            "scipy",
            "pandas",
            "feather-format",
            "biopython",
            "pysam",
            "tqdm",
            "matplotlib-scalebar"
            ],
        entry_points={
            "console_scripts": [
                "hic_hiker = hic_hiker.main:main"
                ]
            }
        )
