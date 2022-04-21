#!/usr/bin/env python
import glob
import setuptools

scripts = []
for s in glob.glob("src/scripts/*"):
    scripts.append(s)

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="fbilr",
    version="1.0.2",
    author="Zonggui Chen",
    author_email="ggchenzonggui@qq.com",
    description="Find barcode in long reads",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://gitee.com/chenzonggui/fbilr",
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Unix"
    ],
    scripts=scripts,
    install_requires=['edlib', 'biopython', 'numpy', 'matplotlib', 'pygz'],
    test_suite="tests",
    python_requires=">=3.6",
)
