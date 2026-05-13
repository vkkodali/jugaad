# jugaad

`jugaad` is a Python command line toolkit for small bioinformatics utilities.

## Installation

Install the CLI from a local checkout:

```sh
git clone git@github.com:vkkodali/jugaad.git
cd jugaad
uv tool install .
jugaad --help
```

For editable development, use `uv sync` and run commands with `uv run`.

## Development

Install the project and development tools with uv:

```sh
uv sync
uv run pre-commit install --hook-type pre-commit --hook-type commit-msg
```

Run the CLI:

```sh
uv run jugaad --help
uv run jugaad gff3_to_introns --help
uv run jugaad gff3_to_introns -f gff3 -i input.gff3 -o introns.tsv
uv run jugaad gff3_to_introns -f gtf --splice_structures < input.gtf
```

Run checks:

```sh
uv run pytest
uv run ruff check .
uv run ruff format --check .
uv run pre-commit run --all-files
```

Format code when needed:

```sh
uv run ruff format .
```

## Releases

This repo uses Conventional Commits and `python-semantic-release`.

Every valid commit on `main` produces at least a patch bump. `feat` commits
can bump minor versions, and commits with `BREAKING CHANGE:` can bump major
versions.

The repository starts at `0.0.0` until the first semantic-release run creates
the first release tag.

The release workflow runs tests, linting, formatting, then:

```sh
uv run semantic-release version
```
