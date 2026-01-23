"""
Setup script for PFAScreen package
"""

from setuptools import setup, find_packages

with open("requirements.txt", "r", encoding="utf-8") as f:
    requirements = [line.strip() for line in f if line.strip() and not line.startswith("#")]

setup(
    name="pfascreen",
    version="1.0.0",
    description="PFAS screening and analysis package",
    author="Your Name",
    author_email="your.email@example.com",
    packages=find_packages(),
    python_requires=">=3.7",
    install_requires=requirements,
    include_package_data=True,
    package_data={
        "pfascreen": ["data/*", "*.yaml"],
    },
    entry_points={
        "console_scripts": [
            "pfascreen-demo=pfascreen.development.demo_pfascreen:main",
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
)
