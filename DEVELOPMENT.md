## Local development setup

After cloning the repository, navigate to the repository root directory
and follow the steps below to set up your development environment.

Create and activate a new environment with the dependencies:

```
conda env create -f environment-dev.yml
conda activate qtrack-dev
```

Editable-install the package:

```
pip install -e . --no-deps
```

Set up the automatic code formatting and linting Git hooks:

1. Install [pre-commit](https://pre-commit.com/):
   ```
   conda install -c conda-forge pre-commit
   ```
   or (if not using Conda)
   ```
   pip install pre-commit
   ```

2. Install the hooks:
   ```
   pre-commit install
   ```

## Running tests

```
pytest
```

## Publish a new version

1. Update/check `__version__` in `qtrack/__init__.py`
   as this is used to set the version for the published package.
   The syntax is described [here](https://packaging.python.org/en/latest/specifications/version-specifiers/#version-specifiers).

2. Use `flit` to build and publish the package:
   ```
   flit publish
   ```

3. Upon successful upload to PyPI, tag the current state of the repository.
   The tag should be named `v<version>`, where `<version>` is the version from step 1.
   For example:
   ```
   git tag -a v0.0.1 -m "v0.0.1"
   git push origin v0.0.1
   ```
