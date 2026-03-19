#!/usr/bin/env python3
"""
Convert Manual.md from Pandoc syntax to MyST/Sphinx syntax
"""
import re

def convert_pandoc_to_myst(content):
    """Convert Pandoc-specific syntax to MyST"""

    # Remove YAML frontmatter (Sphinx doesn't need it)
    content = re.sub(r'^---\n.*?\n---\n', '', content, flags=re.DOTALL)

    # Fix image paths: docs/img/ -> img/
    content = re.sub(r'\(docs/img/', '(img/', content)

    # Convert PDF image references to PNG
    content = re.sub(r'\.pdf\)', '.png)', content)

    # Convert figure labels with width: ![caption](path){#fig:name width=X%}
    def replace_figure_with_width(match):
        caption = match.group(1)
        path = match.group(2)
        name = match.group(3)
        width = match.group(4) if match.lastindex >= 4 and match.group(4) else None

        result = f"```{{figure}} {path}\n:name: fig-{name}\n"
        if width:
            result += f":width: {width}\n"
        result += f"\n{caption}\n```"
        return result

    content = re.sub(r'!\[(.*?)\]\((.*?)\)\{#fig:(\w+)(?:\s+width=(\d+%))?\}',
                     replace_figure_with_width, content)

    # Convert figure references: [@Fig:name] -> {numref}`fig-name`
    content = re.sub(r'\[@Fig:(\w+)\]', r'{numref}`fig-\1`', content)

    # Convert section references: [@Sec:name] -> {ref}`sec-name`
    content = re.sub(r'\[@Sec:(\w+)\]', r'{ref}`sec-\1`', content)

    # Convert table references: [@Tbl:name] -> {numref}`tbl-name`
    content = re.sub(r'\[@Tbl:(\w+)\]', r'{numref}`tbl-\1`', content)

    # Convert section labels: # Title {#sec:name} -> (sec-name)=\n# Title
    content = re.sub(r'^(#{1,6})\s+(.*?)\s+\{#sec:(\w+)\}',
                     r'(sec-\3)=\n\1 \2', content, flags=re.MULTILINE)

    # Convert table labels: {#tbl:name} ->  :name: tbl-name
    # This is approximate and may need manual fixing

    # Convert citations: [@author] -> {cite}`author` for sphinxcontrib-bibtex
    # Handle single citations
    content = re.sub(r'\[@(\w+)\]', r'{cite}`\1`', content)

    # Handle multiple consecutive citations: [@author1][@author2] -> {cite}`author1,author2`
    def merge_citations(match):
        # Extract all citation keys
        keys = re.findall(r'\{cite\}`(\w+)`', match.group(0))
        if len(keys) > 1:
            return '{cite}`' + ','.join(keys) + '`'
        return match.group(0)

    content = re.sub(r'(?:\{cite\}`\w+`\s*)+', merge_citations, content)

    return content

if __name__ == '__main__':
    # Read the original manual
    with open('source/manual.md', 'r', encoding='utf-8') as f:
        content = f.read()

    # Convert
    converted = convert_pandoc_to_myst(content)

    # Write back
    with open('source/manual.md', 'w', encoding='utf-8') as f:
        f.write(converted)

    print("Conversion complete!")
    print("Note: Some manual fixes may still be needed for:")
    print("- Citations (converted to simple [Author] format)")
    print("- Complex tables with labels")
    print("- PDF images (may need conversion to PNG)")
