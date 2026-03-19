#!/usr/bin/env python3
"""
Fix remaining Pandoc syntax in manual.md:
1. Convert equation labels {#eq:name} to MyST math directives
2. Convert table labels {#tbl:name} to proper table captions
3. Convert CoolProp citations to proper {cite} format
"""

import re

def fix_equation_labels(content):
    """Convert inline equations with labels to MyST math blocks."""
    # Pattern: $$ equation $$ {#eq:label}
    # All on single line, so no DOTALL flag needed
    # Pattern breakdown:
    # \$\$ - opening $$
    # \s+ - whitespace
    # (.*?) - equation content (non-greedy, single line)
    # \s+ - whitespace
    # \$\$ - closing $$
    # \s* - optional whitespace
    # \{#eq:([a-zA-Z0-9_-]+)\} - label
    pattern = r'\$\$\s+(.*?)\s+\$\$\s*\{#eq:([a-zA-Z0-9_-]+)\}'

    def replace_equation(match):
        equation = match.group(1).strip()
        label = match.group(2)
        # Convert to MyST math directive
        return f"```{{math}}\n:label: eq-{label}\n\n{equation}\n```"

    content = re.sub(pattern, replace_equation, content)
    return content

def fix_table_labels(content):
    """Convert table labels to proper MyST table captions."""
    # Pattern: : Caption text {#tbl:label}
    pattern = r'^: (.*?) \{#tbl:([a-zA-Z0-9_-]+)\}$'

    def replace_table(match):
        caption = match.group(1).strip()
        label = match.group(2)
        # Keep caption but remove Pandoc label syntax
        # MyST will auto-generate labels from table content
        return f": {caption}"

    content = re.sub(pattern, replace_table, content, flags=re.MULTILINE)
    return content

def fix_coolprop_citations(content):
    """Convert CoolProp DOI references to proper citations."""
    # Pattern: [@doi:10.1021/ie4033999] -> {cite}`coolprop`
    content = content.replace('[@doi:10.1021/ie4033999]', '{cite}`coolprop`')
    return content

def main():
    input_file = '/home/anra/github/HydDown/docs/source/manual.md'

    print("Reading manual.md...")
    with open(input_file, 'r', encoding='utf-8') as f:
        content = f.read()

    print("Fixing equation labels...")
    content = fix_equation_labels(content)

    print("Fixing table labels...")
    content = fix_table_labels(content)

    print("Fixing CoolProp citations...")
    content = fix_coolprop_citations(content)

    print("Writing fixed manual.md...")
    with open(input_file, 'w', encoding='utf-8') as f:
        f.write(content)

    print("Done! Fixed:")
    print("  - Equation labels (converted to MyST math directives)")
    print("  - Table labels (removed Pandoc syntax)")
    print("  - CoolProp citations (converted to {cite})")

if __name__ == '__main__':
    main()
