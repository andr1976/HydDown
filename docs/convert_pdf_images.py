#!/usr/bin/env python3
"""
Convert PDF images to PNG format for better web compatibility
"""
import fitz  # PyMuPDF
from pathlib import Path

def convert_pdf_to_png(pdf_path, output_path=None, dpi=300):
    """Convert a PDF file to PNG image"""
    if output_path is None:
        output_path = pdf_path.with_suffix('.png')

    # Open the PDF
    doc = fitz.open(pdf_path)

    # Get the first page (assuming single-page PDFs)
    page = doc[0]

    # Calculate zoom for desired DPI (72 is default)
    zoom = dpi / 72
    mat = fitz.Matrix(zoom, zoom)

    # Render page to pixmap
    pix = page.get_pixmap(matrix=mat)

    # Save as PNG
    pix.save(output_path)

    doc.close()
    print(f"Converted: {pdf_path.name} -> {output_path.name}")

if __name__ == '__main__':
    # Find all PDF images in source/img
    img_dir = Path('source/img')
    pdf_files = list(img_dir.glob('*.pdf'))

    if not pdf_files:
        print("No PDF files found in source/img/")
    else:
        print(f"Found {len(pdf_files)} PDF files to convert:")
        for pdf_file in pdf_files:
            try:
                convert_pdf_to_png(pdf_file)
            except Exception as e:
                print(f"Error converting {pdf_file.name}: {e}")

        print("\nConversion complete!")
