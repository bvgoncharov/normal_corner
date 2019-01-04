import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="normal_corner",
    version="0.0.1",
    author="Boris Goncharov",
    author_email="goncharov.boris@physics.msu.ru",
    description="Analytical corner plot for Gaussian distributions",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/bvgoncharov/normal_corner",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=["numpy", "scipy", "matplotlib"],
)
