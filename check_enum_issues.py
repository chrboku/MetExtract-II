#!/usr/bin/env python3
"""
Search for potential PySide6 enum issues across the codebase
"""

import os
import re
import glob


def find_potential_enum_issues():
    """Find potential enum issues in Python files"""

    # Common Qt methods that require enum values
    enum_patterns = [
        (r"setVerticalScrollBarPolicy\(([0-9]+)\)", "ScrollBarPolicy enum needed"),
        (r"setHorizontalScrollBarPolicy\(([0-9]+)\)", "ScrollBarPolicy enum needed"),
        (r"setLineWrapMode\(([0-9]+)\)", "LineWrapMode enum needed"),
        (r"setAlignment\(([0-9]+)\)", "Alignment enum needed"),
        (r"setCheckState\(([0-2])\)", "CheckState enum needed"),
        (r"setOrientation\(([0-2])\)", "Orientation enum needed"),
        (r"setPolicy\(([0-9]+)\)", "Policy enum needed"),
    ]

    python_files = []
    # Find all Python files in the project
    for root, dirs, files in os.walk("."):
        for file in files:
            if file.endswith(".py"):
                python_files.append(os.path.join(root, file))

    issues_found = []

    for filepath in python_files:
        try:
            with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
                content = f.read()

            for pattern, description in enum_patterns:
                matches = re.finditer(pattern, content)
                for match in matches:
                    line_num = content[: match.start()].count("\n") + 1
                    issues_found.append(
                        {
                            "file": filepath,
                            "line": line_num,
                            "match": match.group(0),
                            "description": description,
                        }
                    )
        except Exception as e:
            print("Error reading %s: %s" % (filepath, e))

    return issues_found


if __name__ == "__main__":
    print("Searching for potential PySide6 enum issues...")
    issues = find_potential_enum_issues()

    if issues:
        print("\nFound %d potential enum issues:" % len(issues))
        for issue in issues:
            print(
                "  %s:%d - %s (%s)"
                % (issue["file"], issue["line"], issue["match"], issue["description"])
            )
    else:
        print("\nNo obvious enum issues found!")

    print("\nSearch complete.")
