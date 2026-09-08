# 0001: Figures Generation

## Status

Accepted

## Context

The lab needs a standardized approach for generating figures using Python scripts. This ensures consistency, reproducibility, and maintainability of visual outputs.

## Decision

We will use matplotlib as the primary library for creating figures in Python scripts.

We will not use seaborn based on analysis of existing scripts and the need for more direct control over plot elements.

All plots will follow a consistent style guide to ensure uniform appearance across different scripts.

Versioned environments will be used to ensure reproducibility of figures.

## Consequences

- All new Python scripts generating figures will use matplotlib.
- Existing scripts using seaborn will be migrated to matplotlib if necessary.
- A style guide document will be created to standardize visual appearance (colors, fonts, sizes, etc.).
- Versioned environments will be enforced for all figure generation scripts.
