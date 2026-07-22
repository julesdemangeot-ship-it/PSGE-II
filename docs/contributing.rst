Contributing Guide
==================

Thank you for your interest in contributing to PSGE-II!

.. contents::
   :local:
   :depth: 2

Development Setup
-----------------

1. **Clone the repository**::

      git clone https://github.com/julesdemangeot-ship-it/PSGE-II.git
      cd PSGE-II

2. **Create a virtual environment** and install the package in editable mode
   with all development dependencies::

      python -m venv .venv
      source .venv/bin/activate          # Linux / macOS
      # .venv\Scripts\activate           # Windows

      pip install -e ".[dev,docs]"

3. **Verify the installation**::

      python -m pytest tests/ -v

Code Style
----------

PSGE-II enforces consistent formatting using:

* **black** (line length 100) — code formatter
* **isort** (black-compatible profile) — import ordering
* **flake8** — linting

Run all formatters before submitting a pull request::

   black src/ tests/
   isort src/ tests/
   flake8 src/ tests/

Type Annotations
~~~~~~~~~~~~~~~~

All public functions must carry complete type annotations.  Use ``mypy`` to
verify::

   mypy src/psge

Docstrings
~~~~~~~~~~

Use **NumPy-style** docstrings for all public symbols.  Example::

   def gram_matrix(points: np.ndarray) -> np.ndarray:
       """Compute Gram matrix from a set of points.

       Parameters
       ----------
       points : np.ndarray
           Array of shape ``(n, d)`` where *n* is the number of points
           and *d* is the spatial dimension.

       Returns
       -------
       np.ndarray
           Gram matrix of shape ``(n, n)`` with
           ``G[i,j] = <p_i - p_0, p_j - p_0>``.
       """

Running Tests
-------------

Run the full test suite with coverage::

   pytest tests/ --cov=psge --cov-report=term-missing

Aim for **≥ 80 %** line coverage on all implemented modules.

Building the Documentation
--------------------------

From the ``docs/`` directory::

   make html          # Linux / macOS
   make.bat html      # Windows

The generated HTML will be in ``docs/_build/html/``.  Open
``docs/_build/html/index.html`` in a browser to preview.

Pull Request Guidelines
-----------------------

1. **Create a feature branch** from ``main``::

      git checkout -b feature/my-improvement

2. **Write tests** for any new functionality (place them in ``tests/``).
3. **Update the documentation** if the public API changes.
4. **Ensure all checks pass** before opening a PR:

   * ``pytest tests/``
   * ``black --check src/ tests/``
   * ``flake8 src/ tests/``
   * ``mypy src/psge``

5. Submit the pull request with a clear description of the changes.

Version History
---------------

* **v1.2-dev** — Intrinsic formulation (Cayley-Menger based) in progress.
* **v1.1** — Stable extrinsic Euclidean engine.
* **v1.0** — Initial release.
