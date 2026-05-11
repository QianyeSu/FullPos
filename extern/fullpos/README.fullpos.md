# OpenIFS FULLPOS Source Snapshot

This directory contains a local snapshot of selected OpenIFS source directories
used as references for the standalone `fullpos` Python package.

Source root:

```text
F:\openifs-main
```

Copied directories:

```text
F:\openifs-main\ifs-source\arpifs\fullpos  -> upstream\arpifs\fullpos
F:\openifs-main\ifs-source\arpifs\interpol -> upstream\arpifs\interpol
F:\openifs-main\ifs-source\arpifs\module   -> upstream\arpifs\module
F:\openifs-main\ifs-source\trans           -> upstream\trans
F:\openifs-main\ifs-source\contrib         -> upstream\contrib
```

The copied upstream files should be treated as a source snapshot. Package-owned
Python, C, and Fortran wrappers should live under `src/fullpos/`, not inside
`upstream/`.
