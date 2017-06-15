# Development
* [Project Structure](#project-structure)
 * [Code Style](#code-style)
   * [Python](#python)
     * [PEP 8](#pep-8)
     * [Line Length](#line-length)
     * [Multi\-line Literals and Function Calls](#multi-line-literals-and-function-calls)

## Project Structure
The main directories within the source directory are:

* **`batches`**: Contains the input and output files from analysis. This is mostly the domain of users and not developers. Refer to the [batches documentation](batches.md) for more information.
* **`designs`**: Contains the design directories, which are also for user modification, except for the built in ALL and CARDIOM designs, which come with Cpipe. In addition, the genelists directory contains shared region information that might need occasional updating. See the [designs documentation](designs.md) for more information.
* **`docs`**: Contains documentation in markdown format.
* **`lib`**: Contains the Python source code for Cpipe, as a python package. Python-specific tests and scripts are also kept here
* **`pipeline`**: Contains the Groovy source code, namely the core of the pipeline, including tests, and the main pipeline entry point, pipeline.groovy
* **`scripts`**: Interpreted scripts written in languages other than python, that are used by the pipeline or as utility scripts. Includes bash scripts, groovy scripts etc, with the relevant shebang line. Note that these scripts should **not** have a file extension, because this means the binary created from these will also have that extension, which does not meet the unix standard practise. These scripts are copied into the local bin directory by the install process.
* **`tasks`**: Contains tasks written in Python using the [doit](http://pydoit.org/contents.html) automation framework. Currently only used for installation.
* **`tmpdata`**: Used instead of /tmp for the install process, as there is likely more storage in the drive Cpipe is being installed into, rather than the root filesystem.
* **`tools`**: Will contain the tools installed during the installation process.
* **`data`**: Will contain the genomic reference data installed during the installation process.


## Code Style
### Python
#### PEP 8
By default Cpipe uses the [PEP 8](https://www.python.org/dev/peps/pep-0008/) style guide for all python code, except
where the two contradict, in which case the Cpipe style guide takes priority.

#### Line Length
Unlike PEP 8, we prefer a 120 character limit per line. This is the maximum that can be displayed without overflow on a 
github review, and is the default line length in JetBrains IDEs.

#### Multi-line Literals and Function Calls
For any function call, list literal, dictionary literal, or similar comma separated construct, if all arguments/entries
cannot fit on the same line, put each entry on an individual line with a single extra indentation, and then put the 
closing bracket/brace on a new line at the same level of indentation as the opening line.

For instance, this dictionary literal is bad because it is too long for one single line: 
```python
results.append({'interval': candidate, 'distance': 0, 'coding_intersect': coding_intersect, 'exon_codon_positions': exon_codon_positions, 'rank': exon_rank})
```

Instead, it should be written as:
```python
results.append({
    'interval': candidate,
    'distance': 0,
    'coding_intersect': coding_intersect,
    'exon_codon_positions': exon_codon_positions,
    'rank': exon_rank
})
```

Note that the closing section (`})`) is on the same level of indentation as the opening section (`results.append({`), 
and that each entry has one extra level of indentation (another 4 spaces).
