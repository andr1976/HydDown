# HydDown Documentation

This directory contains the Sphinx documentation for HydDown.

## Building the Documentation Locally

### Prerequisites

Install the required dependencies:

```bash
pip install -r ../requirements.txt
```

### Build HTML Documentation

```bash
cd docs
make html
```

The built documentation will be in `build/html/`. Open `build/html/index.html` in your browser to view it.

### Other Build Formats

```bash
make latexpdf  # Build PDF documentation
make epub      # Build EPUB documentation
make clean     # Clean build directory
```

## Documentation Structure

- `source/` - Source files for documentation
  - `index.rst` - Main documentation index
  - `manual.md` - Migrated from Manual.md (comprehensive user guide)
  - `installation.rst` - Installation instructions
  - `quickstart.rst` - Quick start guide
  - `examples.rst` - Usage examples
  - `citing.rst` - Citation information
  - `changelog.rst` - Version history
  - `api/` - API reference documentation
    - `hdclass.rst` - Main calculation engine
    - `transport.rst` - Heat and mass transfer
    - `fire.rst` - Fire heat load modeling
    - `thermesh.rst` - 1-D heat conduction
    - `validator.rst` - Input validation
    - `materials.rst` - Material properties
  - `conf.py` - Sphinx configuration

## Automatic Deployment

Documentation is automatically built and deployed to GitHub Pages when changes are pushed to the main branch. The workflow is defined in `.github/workflows/docs.yml`.

View the live documentation at: https://andr1976.github.io/HydDown/

## Contributing to Documentation

To contribute to the documentation:

1. Edit the appropriate `.rst` or `.md` files in `source/`
2. Build locally to test: `make html`
3. Commit and push changes
4. Documentation will be automatically rebuilt and deployed

## Documentation Features

- **Read the Docs Theme**: Modern, responsive theme
- **MyST Parser**: Supports both Markdown and reStructuredText
- **Autodoc**: Automatic API documentation from docstrings
- **Napoleon**: Support for Google and NumPy style docstrings
- **MathJax**: LaTeX math rendering
- **Intersphinx**: Links to Python, NumPy, SciPy, and Matplotlib documentation
