# MetExtract II - Python 3 Version

A comprehensive metabolite extraction and analysis tool, migrated from Python 2.7 to Python 3 with modern dependency management using UV.

## Migration Notes

This version has been migrated from Python 2.7 to Python 3.8+ with the following major changes:

### Key Changes Made:
- **Python 3 Compatibility**: All Python 2.7 syntax updated to Python 3
- **PyQt4 → PyQt5**: GUI framework updated for modern compatibility
- **Print Statements**: All `print` statements converted to `print()` function calls
- **Dictionary Methods**: `.has_key()` replaced with `in` operator, `.iteritems()` with `.items()`
- **String Handling**: Unicode handling updated for Python 3
- **Integer Division**: Updated `/` to `//` where floor division intended
- **Exception Handling**: Modern `except Exception as e:` syntax
- **Import System**: Relative imports updated to absolute imports where needed

### New Features:
- **UV Package Management**: Modern, fast dependency resolution and virtual environment management
- **Type Hints**: Gradual addition of type annotations
- **Modern Setup**: Uses `pyproject.toml` instead of `setup.py`
- **Code Quality**: Black formatting and linting support

## Installation

### Prerequisites
- Python 3.8 or higher
- UV package manager

### Install UV
```bash
# On Windows
powershell -c "irm https://astral.sh/uv/install.ps1 | iex"

# On macOS/Linux
curl -LsSf https://astral.sh/uv/install.sh | sh
```

### Install MetExtract II
```bash
# Clone the repository
git clone <repository-url>
cd PyMetExtract3

# Create virtual environment and install dependencies
uv venv
uv pip sync requirements.txt

# Or install in development mode
uv pip install -e .
```

### Activate Environment
```bash
# On Windows
.venv\Scripts\activate

# On macOS/Linux
source .venv/bin/activate
```

## Usage

```bash
# Main application
uv run metextract

# Individual modules
uv run mexract
uv run fragextract
uv run fticrmodule
```

## Development

### Code Formatting
```bash
uv run black .
```

### Linting
```bash
uv run flake8 .
```

### Testing
```bash
uv run pytest
```

## Migration Status

✅ **Completed:**
- Python 3 syntax conversion
- PyQt4 → PyQt5 migration
- UV package management setup
- Print statement conversion
- Dictionary method updates
- Basic error handling updates

🔄 **In Progress:**
- Type hint additions
- Test suite updates
- Documentation updates

📋 **TODO:**
- Performance optimization
- Modern Python idiom adoption
- Enhanced error handling
- CI/CD pipeline setup

## Original License

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.