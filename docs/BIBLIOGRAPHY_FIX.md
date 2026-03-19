# Bibliography Fix - Added School Field

## Changes Made

Fixed two thesis entries in `docs/source/references.bib` to use the correct BibTeX field `school` instead of `institution`.

### 1. Fixed iskov Entry (Master's Thesis)

**Before:**
```bibtex
@mastersthesis{iskov,
  ...
  institution = {Aalborg University},  ❌ Wrong field
  ...
}
```

**After:**
```bibtex
@mastersthesis{iskov,
  ...
  school = {Aalborg University},  ✅ Correct field
  ...
}
```

**Result in Bibliography:**
```
Jacob Gram Iskov Eriksen and Michael Skov Bjerre. Analysis of the
application and sizing of pressure safety valves for fire protection
on offshore oil and gas installations. Master's thesis, Aalborg
University, 2015.
```

### 2. Fixed WONG Entry (PhD Thesis)

**Before:**
```bibtex
@phdthesis{WONG,
  ...
  institution = {Department of Chemical Engineering,
                 University College London},  ❌ Wrong field
  ...
}
```

**After:**
```bibtex
@phdthesis{WONG,
  ...
  school = {University College London},  ✅ Correct field
  ...
}
```

**Result in Bibliography:**
```
Shan Meng Angela Wong. DEVELOPMENT OF A MATHEMATICAL MODEL FOR
BLOWDOWN OF VESSELS CONTAINING MULTICOMPONENT HYDROCARBON MIXTURES.
PhD thesis, University College London, 1998.
```

## Why This Matters

### BibTeX Field Requirements

For thesis entries in BibTeX:
- `@mastersthesis` and `@phdthesis` require the `school` field (not `institution`)
- `@techreport` uses `institution` field
- Using wrong field name can cause:
  - Bibliography formatting issues
  - Missing information in citations
  - Build warnings or errors

### Proper Field Usage

| Entry Type | Required Field | Description |
|------------|----------------|-------------|
| `@phdthesis` | `school` | University granting the PhD |
| `@mastersthesis` | `school` | University granting the Master's |
| `@techreport` | `institution` | Organization issuing the report |
| `@book` | `publisher` | Publishing company |

## Build Status

After fix:
- ✅ Build succeeded with 58 warnings (down from 60)
- ✅ Both thesis entries display correctly
- ✅ School names properly shown in bibliography

## Testing

Verified that both entries now appear correctly:

```bash
cd docs
make clean
make html
# Check build/html/manual.html for correct formatting
```

## All Thesis Entries in Bibliography

After this fix, all thesis entries in the bibliography are properly formatted:

1. **iskov** - Master's thesis, Aalborg University ✅
2. **WONG** - PhD thesis, University College London ✅

Both use the correct `school` field and display properly in the bibliography.

## Summary

**What was fixed**: Changed `institution` → `school` for thesis entries
**Files modified**: `docs/source/references.bib`
**Entries fixed**: `iskov`, `WONG`
**Result**: Proper BibTeX formatting, correct bibliography display

Ready to commit with the GitHub Actions fix!
