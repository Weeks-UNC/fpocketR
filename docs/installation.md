# Installation

## Recommended: Conda + pip

1. **Install Conda**
   If you don’t have conda, follow the [official installation guide](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).

2. **Create and activate a new environment with fpocket and Python 3.11**
   ```bash
   conda create -n fpocketR -c conda-forge fpocket=4.0.3 python=3.11
   conda activate fpocketR
   ```

3. **Install fpocketR and dependencies from PyPI**
   ```bash
   pip install fpocketR
   ```
   This will install all required Python dependencies.
   **Note:** The `fpocket` binary is installed via conda, not pip.

---

## Testing your installation

After installing, you can verify your setup by running the test suite:

1. Install the testing tools:
   ```bash
   pip install 'fpocketR[test]'
   ```
2. Find the fpocketR install location:
   ```bash
   python -c "import fpocketR; print(fpocketR.__file__)"
   ```
3. Run `pytest` in the fpocketR source directory:
   ```bash
   pytest /path/to/fpocketR/
   # or, if you are in the source directory:
   pytest
   ```

If all tests pass, your installation is working correctly.

---

## Alternative: Conda Constructor Installer

A one-step installer can be provided using [conda constructor](https://github.com/conda/constructor).
(Instructions and download link will be added here when available.)

---

**Notes:**
- For Windows users, use WSL (Windows Subsystem for Linux) for best compatibility. [Guide to installing WSL and Ubuntu](https://www.freecodecamp.org/news/how-to-install-wsl2-windows-subsystem-for-linux-2-on-windows-10/)
- For MacOS users: fpocketR is not compatible with arm-based M1/M2 processors (only Intel/x86).
