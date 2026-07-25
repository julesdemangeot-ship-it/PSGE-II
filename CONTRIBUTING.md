# Contributing to PSGE-II

Thank you for your interest in contributing to PSGE-II!

## Table of Contents

- [Reporting Issues](#reporting-issues)
- [Submitting a Pull Request](#submitting-a-pull-request)
- [Closing Issues via Pull Requests](#closing-issues-via-pull-requests)
- [Development Setup](#development-setup)
- [Running Tests](#running-tests)
- [Code Style](#code-style)
- [Branching Strategy](#branching-strategy)

---

## Reporting Issues

Before opening a new issue, please search existing issues to avoid duplicates.

When filing a bug report, include:
- Python version and OS
- Minimal reproducible example
- Expected vs. actual behavior
- Full traceback if applicable

---

## Submitting a Pull Request

1. Fork the repository and create a feature branch from `main` (or `develop/v1.2` for v1.2 work).
2. Make your changes with clear, focused commits.
3. Ensure all tests pass and coverage stays above 80%:
   ```bash
   pytest tests/ --cov=psge --cov-fail-under=80
   ```
4. Open a PR targeting the appropriate base branch.
5. Fill in the PR description — see the section below on **closing issues**.

---

## Closing Issues via Pull Requests

Use one of GitHub's [closing keywords](https://docs.github.com/en/issues/tracking-your-work-with-issues/linking-a-pull-request-to-an-issue)
in the PR **description** (not in a commit message) to automatically close the
linked issue when the PR is merged:

```
Closes #<issue-number>
```

Other accepted keywords: `Fixes #N`, `Resolves #N`.

**Example PR description:**

```
## Summary
Add Cayley-Menger determinant guard against near-degenerate simplices.

Closes #3
```

**Rules:**
- One `Closes #N` per issue resolved — list them on separate lines for multiple issues.
- Place the keyword in the PR body, not in a commit message, so GitHub picks it up reliably.
- A single PR may close several issues (e.g. `Closes #3\nCloses #7`).
- If the work is partial, use `Related to #N` instead to avoid premature closure.

---

## Development Setup

```bash
git clone https://github.com/julesdemangeot-ship-it/PSGE-II.git
cd PSGE-II
pip install -e ".[dev]"
```

---

## Running Tests

```bash
# Full test suite with coverage
pytest tests/ --cov=psge --cov-report=term-missing --cov-fail-under=80 -v

# Quick run
pytest tests/
```

---

## Code Style

- Follow [PEP 8](https://peps.python.org/pep-0008/).
- Docstrings use the NumPy convention.
- Keep numerical tolerances dimensionally consistent (prefer relative/adimensional thresholds over hard-coded absolute values).

---

## Branching Strategy

| Branch | Purpose |
|--------|---------|
| `main` | Integration branch — PRs must pass CI |
| `stable/v1.1` | Protected stable reference — no direct pushes |
| `develop/v1.2` | Active feature development |
| `feature/*` | Short-lived feature branches |
| `research/*` | Exploratory / long-horizon work |

Target `develop/v1.2` for v1.2 features; target `main` for documentation, CI, and tooling changes.
