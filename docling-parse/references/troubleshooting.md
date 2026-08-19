# Docling Troubleshooting Guide

Common failure modes, their causes, and fixes for fleet use.

---

## Installation Issues

### `ImportError: No module named 'docling'`
**Cause:** Docling not installed in the active Python environment.
**Fix:**
```bash
pip install docling
# Verify:
python -c "import docling; print(docling.__version__)"
```

### `Python 3.9 not supported`
**Cause:** Docling 2.70.0+ requires Python 3.10+.
**Fix:**
```bash
python --version   # Check current version
# If 3.9, use pyenv or conda to switch:
pyenv install 3.11.0
pyenv local 3.11.0
pip install docling
```

### Slow first install (>10 minutes)
**Cause:** Docling has many dependencies (torch, transformers, etc.).
**Fix:** Normal. Use `pip install docling --progress` to see progress. Consider pre-installing on fleet nodes rather than per-run install.

### `ERROR: Could not build wheels for some packages`
**Cause:** Missing build tools (common on fresh Linux nodes).
**Fix:**
```bash
# Ubuntu/Debian:
sudo apt-get install build-essential python3-dev
pip install docling

# macOS:
xcode-select --install
pip install docling
```

---

## Parsing Issues

### Blank or empty output
**Possible causes and fixes:**
1. **Scanned PDF without OCR** — Enable OCR:
   ```bash
   pip install docling[ocr]
   ```
   ```python
   result = parse_pdf("scan.pdf", use_ocr=True)
   ```
2. **Password-protected PDF** — Decrypt first:
   ```bash
   pip install pikepdf
   ```
   ```python
   import pikepdf
   with pikepdf.open("protected.pdf", password="secret") as pdf:
       pdf.save("decrypted.pdf")
   result = parse_pdf("decrypted.pdf")
   ```
3. **Corrupted PDF** — Validate with:
   ```bash
   python -c "import pikepdf; pikepdf.open('file.pdf')"
   ```

### Tables missing from output
**Cause:** Table detection confidence below threshold, or table is an image (scanned).
**Fix:**
```python
# Lower detection threshold
result = parse_pdf("paper.pdf", table_detection_threshold=0.3)

# For scanned tables, enable OCR
result = parse_pdf("paper.pdf", use_ocr=True)
```

### Wrong reading order (columns jumbled)
**Cause:** Complex multi-column layout where layout model confidence is low.
**Fix:**
- Usually improves with newer Docling versions (`pip install --upgrade docling`)
- For critical documents, check with `validate_output.py` and manually verify section order
- Flag affected tables with `[VERIFY]` per First Law

### Wrong column order in tables
**Cause:** Spanning cells (colspan/rowspan) confuse the table structure model. Known Docling limitation.
**Fix:**
- Export to JSON and inspect raw `table_cells` for `col_span` values
- Post-process with pandas if needed
- Flag with `[VERIFY]` in any downstream use

### Unicode/encoding errors in output
**Cause:** Source PDF has non-standard encoding or ligatures.
**Fix:**
```python
# Script handles this automatically; if you see garbled text:
text = text.encode('utf-8', errors='replace').decode('utf-8')
# Or use ftfy library:
pip install ftfy
import ftfy
text = ftfy.fix_text(text)
```

### Title not extracted
**Cause:** PDF metadata doesn't contain title, or first-page layout doesn't have a detectable title.
**Fix:**
- Title falls back to filename (e.g., `paper.pdf` → `"paper"`)
- Set manually: pass `title` kwarg in scripts, or update after parsing

### Figure extraction fails
**Cause:** Write permission issue or `figures/` dir not created.
**Fix:** The `parse_pdf.py` script creates `figures/` automatically. Check:
```bash
ls -la figures/   # Should exist after run
# If permission denied:
chmod 755 .
```

---

## Performance Issues

### Very slow on large PDFs (>50 pages)
**Cause:** Full layout analysis is O(n pages).
**Options:**
1. **Fast mode** (text only, no layout analysis):
   ```python
   result = parse_pdf("big.pdf", fast_mode=True)
   ```
2. **Page range** (not directly in script, use converter options):
   ```python
   from docling.datamodel.document import ConversionResult
   # Limit pages via PdfPipelineOptions — see docling docs
   ```
3. **Split PDF first**:
   ```bash
   pip install pikepdf
   python -c "
   import pikepdf
   pdf = pikepdf.open('big.pdf')
   # Save pages 0-49
   new_pdf = pikepdf.Pdf.new()
   for page in pdf.pages[:50]:
       new_pdf.pages.append(page)
   new_pdf.save('big_part1.pdf')
   "
   ```

### Memory error on batch processing
**Cause:** Many large PDFs loaded into memory simultaneously.
**Fix:**
```bash
# Use smaller batches in parse_batch.py:
python scripts/parse_batch.py --input-dir ./pdfs/ --output-dir ./out/ --batch-size 5
```

### Slow first run (model download)
**Cause:** Docling downloads the Heron layout model and TableFormer on first call (~200–400MB).
**Fix:** Normal behavior. Pre-download on fleet nodes:
```bash
# Pre-download all models:
python -c "from docling.document_converter import DocumentConverter; DocumentConverter()"
# Or via CLI:
docling --download-models
```

---

## OCR-Specific Issues

### `EasyOCR not found` / `pytesseract not found`
**Cause:** OCR extras not installed.
**Fix:**
```bash
pip install docling[ocr]
# For tesseract backend, also install system package:
# macOS:
brew install tesseract
# Ubuntu:
sudo apt-get install tesseract-ocr
```

### OCR produces garbled text
**Cause:** Image quality too low, or wrong language model.
**Fix:**
```python
from docling.datamodel.pipeline_options import PdfPipelineOptions, EasyOcrOptions

ocr_options = EasyOcrOptions(lang=["en"])   # Set correct language
pipeline_options = PdfPipelineOptions(do_ocr=True, ocr_options=ocr_options)
```

### OCR very slow on CPU
**Cause:** OCR models are compute-intensive without GPU.
**Fix:** Route scanned PDF jobs to Buster or Sparq:
```bash
# On Buster (3080 Ti):
python ~/.openclaw/workspace/skills/docling-parse/scripts/parse_pdf.py --input /path/to/scan.pdf --output out.md --format markdown
```

---

## Output Quality Issues

### Output fails `validate_output.py` checks

**Low word count (<100 words):**
- PDF may be mostly figures/tables with little text
- Check if it's a scanned PDF requiring OCR
- Inspect raw PDF: `pdftotext paper.pdf -` (from poppler)

**No section headings detected:**
- Document may use font size/bold rather than structural headings
- Check with `fast_mode=False` (default) to ensure layout analysis ran
- Some PDFs use paragraph-style headings that Docling can't detect — [VERIFY] and proceed

**Tables reported as 0:**
- Enable table detection: ensure `fast_mode=False`
- Check if tables are images (scanned) → use OCR
- Lower `table_detection_threshold` to 0.3

---

## Version/Compatibility Issues

### `AttributeError: 'DoclingDocument' has no attribute 'export_to_dict'`
**Cause:** Very old Docling version.
**Fix:**
```bash
pip install --upgrade docling
```

### API changed between versions
Docling is under active development. Breaking changes can occur between minor versions.
**Fix:**
```bash
pip show docling   # Check current version
# Pin to a known-good version if needed:
pip install "docling==2.x.x"
```

---

## Fleet-Specific Issues

### Different results on different fleet nodes
**Cause:** Different Docling versions installed on different nodes.
**Fix:** Standardize version across fleet:
```bash
# Check version on each node:
for node in host-a host-b host-c; do
  ssh zevo@$node "pip show docling | grep Version"
done
# Install same version everywhere
```

### Model cache not found on node
**Cause:** Models cached for different user, or cache cleared.
**Fix:**
```bash
# Check cache location:
python -c "from docling.utils.model_utils import get_cache_dir; print(get_cache_dir())"
# Re-download:
python -c "from docling.document_converter import DocumentConverter; DocumentConverter()"
```

---

## Reporting Issues

If you encounter a bug not covered here:
1. Record the exact Docling version: `pip show docling`
2. Save the failing PDF (if not sensitive)
3. Note the exact error message and traceback
4. Check [Docling GitHub Issues](https://github.com/docling-project/docling/issues) before opening a new one
5. Add to this troubleshooting guide with Fix once resolved

---

*Reference compiled by Zevo 🦓 | 2026-05-27*
