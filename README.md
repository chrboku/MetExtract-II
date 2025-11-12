# MetExtract II - Python 3 Version

A comprehensive metabolite extraction and analysis tool, migrated from Python 2.7 to Python 3 with modern dependency management using UV.

## Installation

### Optional: Install uv
Please install uv. For instruction see [https://docs.astral.sh/uv/getting-started/installation/](https://docs.astral.sh/uv/getting-started/installation/)

### Clone MetExtract II
```bash
# Clone the repository
git clone https://github.com/chrboku/MetExtract-II
cd MetExtract-II
```

### Setup virtual environment
```bash
uv sync
```

## Quick Start
```bash
uv run python -m src.MExtract
```

## Usage

### Running the Applications

```bash
# Main application (MExtract) - Primary metabolite extraction interface
uv run python -m src.MExtract

# MetExtract II Main Interface - Alternative main interface
uv run python -m src.MetExtractII_Main
```

### Module Descriptions

- **MExtract**: The primary GUI application for metabolite extraction and analysis
- **MetExtractII_Main**: Main interface for MetExtract II functionality
- **FragExtract**: Specialized tool for extracting and analyzing mass spectrometry fragments
- **FTICRModule**: Module for handling Fourier Transform Ion Cyclotron Resonance mass spectrometry data


### Code Formatting
```bash
# Format code with Black
uvx ruff format
```


## Migration Notes

This version has been migrated from Python 2.7 (published version) to Python 3.11+ with the following major changes:

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
