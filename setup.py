from setuptools import setup

setup(
        name="hic_hiker",
        version="0.0.1",
        description="",
        long_description=open('README.md').read(),
        long_description_content_type="text/markdown",
        url="https://github.com/ryought/HiC-Hiker",
        author="ryought",
        author_email="ryonakabayashi@gmail.com",
        license="MIT",
        packages=[
            "hic_hiker",
            ],
        entry_points={
            "console_scripts": [
                "hic_hiker = hic_hiker.main:main"
                ]
            }
        )
