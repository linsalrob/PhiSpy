# PhiSpy Usage in Pyodide

This document provides examples of using PhiSpy in a Pyodide environment (web browser or Node.js).

## Installation in Pyodide

Once PhiSpy is available in the Pyodide distribution or via the pyodide-recipes build:

```python
import micropip
await micropip.install('phispy')
```

## Important Notes

### Dependencies

PhiSpy has the following dependencies:

1. **Available in Pyodide by default:**
   - `biopython` - Built into Pyodide
   - `numpy` - Built into Pyodide
   - `scikit-learn` - Built into Pyodide

2. **Pure Python packages (included in recipe):**
   - `bcbio-gff` - Automatically installed with PhiSpy

3. **Conditional dependency:**
   - `importlib-resources` - Only needed for Python < 3.9 (Pyodide uses Python 3.11+, so not needed)

### C++ Extension

PhiSpy includes a compiled C++ extension (`PhiSpyRepeatFinder`) for efficient repeat sequence finding. This extension is compiled to WebAssembly and works transparently in the browser environment.

## Basic Usage Example

```python
# In a Jupyter notebook or Python REPL in the browser
import micropip
await micropip.install('phispy')

# Import PhiSpy modules
from PhiSpyModules import main
from PhiSpyModules import helper_functions

# The main entry point expects command-line arguments
# In a browser context, you would typically work with the modules directly
# rather than using the command-line interface

print("PhiSpy is ready to use!")
```

## Working with GenBank Files

Since PhiSpy analyzes GenBank files, you'll need to provide genomic data. In a browser environment, you can:

### Option 1: Upload Files

```python
from js import document
from pyodide.ffi import to_js

# Create a file input element
file_input = document.createElement('input')
file_input.type = 'file'
file_input.accept = '.gb,.gbk,.genbank'

async def handle_file(event):
    file = event.target.files.item(0)
    if file:
        # Read file content
        content = await file.text()
        # Save to virtual filesystem
        with open('genome.gb', 'w') as f:
            f.write(content)
        print(f"Loaded {file.name}")

file_input.onchange = handle_file
document.body.appendChild(file_input)
```

### Option 2: Fetch Remote Files

```python
import pyodide.http
import os

# Fetch a test GenBank file
url = "https://raw.githubusercontent.com/linsalrob/PhiSpy/master/tests/Streptococcus_pyogenes_M1_GAS.gb"
response = await pyodide.http.pyfetch(url)
content = await response.string()

# Write to virtual filesystem
with open('test_genome.gb', 'w') as f:
    f.write(content)

print("Downloaded test genome file")
```

## Running PhiSpy Analysis

```python
from Bio import SeqIO
from PhiSpyModules import classification
from PhiSpyModules import helper_functions

# Read the GenBank file
genome_file = 'test_genome.gb'
records = list(SeqIO.parse(genome_file, 'genbank'))

print(f"Loaded {len(records)} record(s)")

# PhiSpy analysis would typically use the main module
# For detailed usage, refer to the PhiSpy documentation
# Note: Full CLI functionality may need adaptation for browser use
```

## Limitations in Browser Environment

When using PhiSpy in a browser via Pyodide, be aware of:

1. **File I/O**: All file operations use Pyodide's virtual filesystem
2. **Memory**: Large genomes may exceed browser memory limits
3. **Performance**: WebAssembly is fast but may be slower than native code
4. **CLI Tools**: Command-line interface features may need adaptation
5. **No Direct Filesystem Access**: Must upload or fetch files explicitly

## Advanced: Creating a Web Interface

```html
<!DOCTYPE html>
<html>
<head>
    <script src="https://cdn.jsdelivr.net/pyodide/v0.29.2/full/pyodide.js"></script>
</head>
<body>
    <h1>PhiSpy Web Interface</h1>
    <input type="file" id="gbFile" accept=".gb,.gbk">
    <button id="analyze">Analyze</button>
    <pre id="output"></pre>
    
    <script type="text/javascript">
        async function main() {
            let pyodide = await loadPyodide();
            await pyodide.loadPackage('micropip');
            await pyodide.runPythonAsync(`
                import micropip
                await micropip.install('phispy')
                print('PhiSpy loaded successfully!')
            `);
            
            document.getElementById('analyze').addEventListener('click', async () => {
                const fileInput = document.getElementById('gbFile');
                const file = fileInput.files[0];
                if (file) {
                    const content = await file.text();
                    pyodide.FS.writeFile('input.gb', content);
                    
                    const result = await pyodide.runPythonAsync(`
                        # Run PhiSpy analysis here
                        "Analysis complete!"
                    `);
                    
                    document.getElementById('output').textContent = result;
                }
            });
        }
        main();
    </script>
</body>
</html>
```

## Performance Tips

1. **Process smaller genomes first** - Test with smaller files before large genomes
2. **Use worker threads** - Consider using Web Workers for heavy computations
3. **Monitor memory** - Large training sets may require more memory
4. **Cache results** - Store intermediate results in browser storage if needed

## Troubleshooting

### Memory Issues with Large Genomes

- Use smaller genome files for testing
- Consider server-side processing for large-scale analysis
- Clear unused data from memory

### Performance is Slow

- This is expected - WebAssembly is generally slower than native code
- For production use, consider running PhiSpy server-side
- The browser version is best for demos, education, and small analyses

## Resources

- [PhiSpy Documentation](https://github.com/linsalrob/PhiSpy)
- [Pyodide Documentation](https://pyodide.org/)
- [BioPython Documentation](https://biopython.org/)

## Citation

If you use PhiSpy in your research, please cite:
```
Akhter, S., Aziz, R. K., & Edwards, R. A. (2012). 
PhiSpy: a novel algorithm for finding prophages in bacterial genomes.
Nucleic Acids Research, 40(16), e126.
```
