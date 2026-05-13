# AI Instructions

- Add new commands as argparse subcommands registered from `jugaad.cli`.
- Keep command modules import-safe: no work should run at import time.
- Prefer small, testable functions for parsing, I/O, and formatting behavior.
