# Cpipe Style Guide
## Python
### PEP 8
By default Cpipe uses the [PEP 8](https://www.python.org/dev/peps/pep-0008/) style guide for all python code, except
where the two contradict, in which case the Cpipe style guide takes priority.

### Line Length
Unlike PEP 8, we prefer a 120 character limit per line. This is the maximum that can be displayed without overflow on a 
github review, and is the default line length in JetBrains IDEs.

### Multi-line Literals and Function Calls
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